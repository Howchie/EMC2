// The chain's side of an exchange with its worker pool (R/chain_pool.R):
// .emc_wpool_subject_ll_queue(), .emc_wpool_iter() over the shared subject
// queue and a round of per-worker requests. Messages are serialized and
// framed here, the shared broadcast is a file the workers read, and the
// completion records and replies are read from the pool's connections. R is
// called only when the pool fails: to warn, to stop stuck workers and to
// compute what no worker answered for.
#include "pool_exchange.h"
#include "pool_transport.h"
#include "iteration_state.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#ifndef _WIN32
#include <dirent.h>
#include <signal.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

using namespace Rcpp;
using emc_iter::flag;
using emc_iter::lget;
using emc_iter::lset;
using emc_iter::package_function;
using emc_iter::set_var;
using emc_iter::var;

void emc_claim_reset(SEXP map, int key);

namespace {

typedef std::chrono::steady_clock Steady;

double seconds_since(const Steady::time_point& t0) {
  return std::chrono::duration<double>(Steady::now() - t0).count();
}

// --- the pool's session state (.emc_pool_state) ---------------------------

SEXP pool_state() { return emc_iter::package_value(".emc_pool_state"); }

void state_add_int(const char* name, int by) {
  SEXP e = pool_state();
  SEXP cur = var(e, name);
  set_var(e, name, Rf_ScalarInteger((cur == R_NilValue ? 0 : Rf_asInteger(cur)) + by));
}

void state_add_real(const char* name, double by) {
  SEXP e = pool_state();
  SEXP cur = var(e, name);
  set_var(e, name, Rf_ScalarReal((cur == R_NilValue ? 0.0 : Rf_asReal(cur)) + by));
}

// .emc_wpool_next_request().
int next_request() {
  SEXP e = pool_state();
  SEXP cur = var(e, "request_seq");
  int n = 1;
  if (cur != R_NilValue) {
    const int last = Rf_asInteger(cur);
    n = last >= INT_MAX - 1 ? 1 : last + 1;
  }
  set_var(e, "request_seq", Rf_ScalarInteger(n));
  return n;
}

// .emc_wpool_deadline().
double deadline() {
  SEXP opt = Rf_GetOption1(Rf_install("emc2.worker_timeout"));
  if (opt != R_NilValue && Rf_length(opt) >= 1) {
    const double v = Rf_asReal(opt);
    if (!ISNAN(v)) return v;
  }
  SEXP h = var(pool_state(), "iter_receive_times");
  if (TYPEOF(h) != REALSXP || XLENGTH(h) == 0) return 600.0;
  std::vector<double> x(REAL(h), REAL(h) + XLENGTH(h));
  std::sort(x.begin(), x.end());
  const size_t m = x.size();
  const double med = (m % 2) ? x[m / 2] : (x[m / 2 - 1] + x[m / 2]) / 2;
  return std::max(120.0, 20 * med);
}

// .emc_wpool_record_receive().
void record_receive(double seconds) {
  SEXP e = pool_state();
  SEXP h = var(e, "iter_receive_times");
  std::vector<double> x;
  if (TYPEOF(h) == REALSXP) x.assign(REAL(h), REAL(h) + XLENGTH(h));
  x.push_back(seconds);
  const size_t from = x.size() > 32 ? x.size() - 32 : 0;
  NumericVector out(x.begin() + from, x.end());
  set_var(e, "iter_receive_times", out);
}

void transport_fail(const std::string& msg) {
  set_var(pool_state(), "last_transport_error", Rf_mkString(msg.c_str()));
}

void degraded(const std::string& msg) {
  static Function warn = package_function(".emc_wpool_degraded");
  warn(msg);
}

// --- serialization ---------------------------------------------------------

void out_char(R_outpstream_t stream, int c) {
  static_cast<std::vector<unsigned char>*>(stream->data)->push_back(
      static_cast<unsigned char>(c));
}

void out_bytes(R_outpstream_t stream, void* buf, int len) {
  std::vector<unsigned char>* out = static_cast<std::vector<unsigned char>*>(stream->data);
  const unsigned char* b = static_cast<const unsigned char*>(buf);
  out->insert(out->end(), b, b + len);
}

// serialize(x, NULL).
void serialize(SEXP x, std::vector<unsigned char>& out) {
  out.clear();
  struct R_outpstream_st stream;
  R_InitOutPStream(&stream, static_cast<R_pstream_data_t>(&out), R_pstream_xdr_format, 3,
                   out_char, out_bytes, NULL, R_NilValue);
  R_Serialize(x, &stream);
}

struct InBytes {
  const unsigned char* p;
  size_t n, at;
};

int in_char(R_inpstream_t stream) {
  InBytes* in = static_cast<InBytes*>(stream->data);
  if (in->at >= in->n) Rf_error("the message ends early");
  return in->p[in->at++];
}

void in_bytes(R_inpstream_t stream, void* buf, int len) {
  InBytes* in = static_cast<InBytes*>(stream->data);
  if (len < 0 || in->at + static_cast<size_t>(len) > in->n) Rf_error("the message ends early");
  std::memcpy(buf, in->p + in->at, len);
  in->at += len;
}

struct Unserialized {
  InBytes in;
  SEXP value;
};

void unserialize_body(void* data) {
  Unserialized* job = static_cast<Unserialized*>(data);
  struct R_inpstream_st stream;
  R_InitInPStream(&stream, static_cast<R_pstream_data_t>(&job->in), R_pstream_any_format,
                  in_char, in_bytes, NULL, R_NilValue);
  job->value = R_Unserialize(&stream);
}

// unserialize(); NULL when the bytes are not a message (R reports why).
SEXP unserialize(const std::vector<unsigned char>& bytes) {
  Unserialized job = {{bytes.data(), bytes.size(), 0}, R_NilValue};
  if (!R_ToplevelExec(unserialize_body, &job)) {
    transport_fail("the reply is not a message");
    return R_NilValue;
  }
  return job.value;
}

// --- the workers -----------------------------------------------------------

// .emc_wpool_pid_alive().
bool pid_alive(int pid) {
#ifdef _WIN32
  return true;
#else
  if (kill(pid, 0) != 0) return false;
  struct stat sb;
  if (stat("/proc", &sb) != 0 || !S_ISDIR(sb.st_mode)) return true;
  char path[64];
  std::snprintf(path, sizeof path, "/proc/%d/stat", pid);
  std::FILE* f = std::fopen(path, "r");
  if (f == NULL) return false;
  char line[1024];
  const bool got = std::fgets(line, sizeof line, f) != NULL;
  std::fclose(f);
  if (!got) return false;
  // After "(comm) ": the state, then the fields up to the flags (the seventh).
  const char* rest = std::strrchr(line, ')');
  if (rest == NULL || rest[1] != ' ') return true;
  rest += 2;
  if (rest[0] == 'Z' || rest[0] == 'X' || rest[0] == 'x') return false;
  const char* at = rest;
  for (int field = 1; field < 7 && at != NULL; ++field) {
    at = std::strchr(at, ' ');
    if (at != NULL) ++at;
  }
  if (at == NULL) return true;
  // A process past exit() still runs while the kernel tears it down.
  const unsigned long long flags = std::strtoull(at, NULL, 10);
  return (flags & 4ULL) == 0;
#endif
}

// A pool's connections and processes for one exchange.
struct Workers {
  int n;
  std::vector<void*> wcs, rcs;
  std::vector<SEXP> jobs;
  std::vector<int> pids;
  std::vector<bool> pid_known;

  explicit Workers(SEXP pool) {
    n = Rf_asInteger(lget(pool, "n"));
    SEXP w = lget(pool, "wcs"), r = lget(pool, "rcs"), j = lget(pool, "jobs");
    wcs.assign(n, NULL);
    rcs.assign(n, NULL);
    jobs.assign(n, R_NilValue);
    pids.assign(n, NA_INTEGER);
    pid_known.assign(n, false);
    for (int k = 0; k < n; ++k) {
      if (TYPEOF(w) == VECSXP && k < XLENGTH(w)) wcs[k] = emc_pool::connection(VECTOR_ELT(w, k));
      if (TYPEOF(r) == VECSXP && k < XLENGTH(r)) rcs[k] = emc_pool::connection(VECTOR_ELT(r, k));
      if (TYPEOF(j) == VECSXP && k < XLENGTH(j)) jobs[k] = VECTOR_ELT(j, k);
    }
  }

  // .emc_wpool_job_pid().
  int pid(int w) {
    if (pid_known[w]) return pids[w];
    pid_known[w] = true;
    SEXP job = jobs[w];
    SEXP own = lget(job, "pid");
    if (own != R_NilValue && Rf_length(own) == 1 && Rf_asInteger(own) != NA_INTEGER) {
      pids[w] = Rf_asInteger(own);
    } else if (job != R_NilValue) {
      static Function job_pid = package_function(".emc_wpool_job_pid");
      pids[w] = Rf_asInteger(job_pid(job));
    }
    return pids[w];
  }

  bool alive(int w) {
    if (jobs[w] == R_NilValue) return true;
    const int id = pid(w);
    return id == NA_INTEGER || pid_alive(id);
  }

  // .emc_wpool_busy() of these workers.
  bool busy(const std::vector<int>& ws) {
    IntegerVector ids(ws.size());
    for (size_t k = 0; k < ws.size(); ++k) ids[k] = jobs[ws[k]] == R_NilValue ? NA_INTEGER : pid(ws[k]);
    static Function is_busy = package_function(".emc_wpool_busy");
    return flag(is_busy(ids));
  }

  // The payload to worker w within `limit` seconds; what went wrong otherwise.
  std::string send(int w, const std::vector<unsigned char>& payload, double limit) {
    const std::string err = emc_pool::send_payload(wcs[w], payload.data(), payload.size(), true,
                                                   limit);
    if (!err.empty()) transport_fail(err);
    return err;
  }

  // Worker w's reply; NULL, with the reason recorded, when there is none.
  SEXP reply(int w, double limit) {
    emc_pool::Watch watch;
    watch.timed = true;
    watch.limit = limit;
    if (jobs[w] != R_NilValue) {
      watch.alive = [this, w]() { return alive(w); };
      watch.busy = [this, w]() { return busy(std::vector<int>(1, w)); };
    }
    std::vector<unsigned char> payload;
    const std::string err = emc_pool::recv_payload(rcs[w], watch, payload);
    if (!err.empty()) {
      transport_fail(err);
      return R_NilValue;
    }
    return unserialize(payload);
  }
};

// .emc_wpool_transport_error().
std::string transport_error(const char* fallback) {
  SEXP msg = var(pool_state(), "last_transport_error");
  if (TYPEOF(msg) == STRSXP && XLENGTH(msg) > 0 && CHAR(STRING_ELT(msg, 0))[0] != '\0')
    return CHAR(STRING_ELT(msg, 0));
  return fallback;
}

// --- the shared broadcast --------------------------------------------------

// A broadcast file, removed when the exchange ends.
struct SharedFile {
  std::string path;
  ~SharedFile() {
    if (!path.empty()) std::remove(path.c_str());
  }
};

// .emc_wpool_shared_write(): the bytes as a new file of the pool's
// directory, complete before it has its name. An empty string on success.
std::string shared_write(SEXP dir, const std::vector<unsigned char>& bytes, std::string& path) {
#ifdef _WIN32
  return "shared broadcast directory is unavailable";
#else
  struct stat sb;
  if (TYPEOF(dir) != STRSXP || XLENGTH(dir) != 1 || stat(CHAR(STRING_ELT(dir, 0)), &sb) != 0 ||
      !S_ISDIR(sb.st_mode))
    return "shared broadcast directory is unavailable";
  if (bytes.size() > static_cast<size_t>(INT_MAX)) return "shared broadcast payload is too large";
  static unsigned long counter = 0;
  char name[96];
  std::snprintf(name, sizeof name, "/shared_%d_%lu.bin", static_cast<int>(getpid()), ++counter);
  const std::string final_path = std::string(CHAR(STRING_ELT(dir, 0))) + name;
  const std::string tmp = final_path + ".tmp";
  std::FILE* f = std::fopen(tmp.c_str(), "wb");
  if (f == NULL) return "cannot open the shared broadcast file";
  const bool wrote = bytes.empty() || std::fwrite(bytes.data(), 1, bytes.size(), f) == bytes.size();
  const bool closed = std::fclose(f) == 0;
  if (!wrote || !closed || std::rename(tmp.c_str(), final_path.c_str()) != 0) {
    std::remove(tmp.c_str());
    return "could not publish shared broadcast file";
  }
  path = final_path;
  return "";
#endif
}

// --- the shared subject queue ---------------------------------------------

struct Claim {
  RObject ptr, path;
};

// The pool's claim counter where the pool serves the shared subject queue.
bool claim_attach(SEXP pool, Claim& claim) {
  SEXP dir = lget(pool, "dir");
  if (!flag(lget(pool, "alive")) || Rf_asInteger(lget(pool, "n")) <= 1 ||
      TYPEOF(dir) != STRSXP || XLENGTH(dir) != 1 || !emc_iter::option_on("emc2.worker_queue"))
    return false;
  const std::string path = std::string(CHAR(STRING_ELT(dir, 0))) + "/claim";
  SEXP cached = var(pool_state(), "claim_map");
  SEXP cached_path = lget(cached, "path");
  RObject map = cached;
  if (TYPEOF(cached_path) != STRSXP || XLENGTH(cached_path) != 1 ||
      path != CHAR(STRING_ELT(cached_path, 0))) {
    static Function attach = package_function(".emc_wpool_claim_attach");
    map = attach(path);
  }
  if (TYPEOF(map) != VECSXP) return false;
  claim.ptr = lget(map, "ptr");
  claim.path = lget(map, "path");
  return TYPEOF(claim.ptr) == EXTPTRSXP;
}

// .emc_wpool_queue_order(): the dispatch order of the subjects in `part`.
// FALSE unless they are the n subjects, each once.
bool queue_order(SEXP part, int n, std::vector<int>& order) {
  if (TYPEOF(part) != VECSXP) return false;
  std::vector<int> subs;
  for (R_xlen_t w = 0; w < XLENGTH(part); ++w) {
    SEXP x = VECTOR_ELT(part, w);
    if (TYPEOF(x) == INTSXP) {
      subs.insert(subs.end(), INTEGER(x), INTEGER(x) + XLENGTH(x));
    } else if (TYPEOF(x) == REALSXP) {
      for (R_xlen_t k = 0; k < XLENGTH(x); ++k) subs.push_back(static_cast<int>(REAL(x)[k]));
    } else if (x != R_NilValue) {
      return false;
    }
  }
  if (static_cast<int>(subs.size()) != n) return false;
  std::vector<bool> held(n + 1, false), placed(n + 1, false);
  for (int s : subs) {
    if (s < 1 || s > n || held[s]) return false;
    held[s] = true;
  }
  order.clear();
  SEXP o = Rf_getAttrib(part, Rf_install("order"));
  if (TYPEOF(o) == INTSXP || TYPEOF(o) == REALSXP) {
    for (R_xlen_t k = 0; k < XLENGTH(o); ++k) {
      const int s = TYPEOF(o) == INTSXP ? INTEGER(o)[k] : static_cast<int>(REAL(o)[k]);
      if (s < 1 || s > n || placed[s]) continue;
      placed[s] = true;
      order.push_back(s);
    }
  }
  for (int s : subs) if (!placed[s]) {
    placed[s] = true;
    order.push_back(s);
  }
  return true;
}

// A reply's subjects, where they are subjects of this round no other reply
// holds.
bool reply_subjects(SEXP subs, int n, const std::vector<bool>& stored, std::vector<int>& out) {
  out.clear();
  if (TYPEOF(subs) != INTSXP && TYPEOF(subs) != REALSXP) return false;
  std::vector<bool> seen(n, false);
  for (R_xlen_t k = 0; k < XLENGTH(subs); ++k) {
    const double v = TYPEOF(subs) == INTSXP
        ? (INTEGER(subs)[k] == NA_INTEGER ? NA_REAL : INTEGER(subs)[k]) : REAL(subs)[k];
    if (ISNAN(v) || v < 1 || v > n) return false;
    const int s = static_cast<int>(v) - 1;
    if (seen[s] || stored[s]) return false;
    seen[s] = true;
    out.push_back(s);
  }
  return true;
}

// a + b of two reject tallies.
SEXP add_counts(SEXP a, SEXP b) {
  if (TYPEOF(a) == INTSXP && TYPEOF(b) == INTSXP && XLENGTH(a) == XLENGTH(b)) {
    Shield<SEXP> out(Rf_duplicate(a));
    for (R_xlen_t k = 0; k < XLENGTH(a); ++k) INTEGER(out)[k] += INTEGER(b)[k];
    return out;
  }
  static Function plus("+");
  return plus(a, b);
}

// --- the completion records ------------------------------------------------

const size_t kRecordBytes = 16;

struct Record {
  int w, generation, request, status;
};

// The pool's completion channel: fixed-size records, read without waiting.
// Bytes of a record that has not fully arrived stay in the channel's buffer.
struct Records {
  void* rc;
  RObject buf;
  unsigned char bytes[kRecordBytes];
  size_t have;

  explicit Records(SEXP done) : rc(emc_pool::connection(lget(done, "rc"))),
                                buf(lget(done, "buf")), have(0) {
    if (TYPEOF(buf) != ENVSXP) return;
    SEXP held = var(buf, "bytes");
    if (TYPEOF(held) == RAWSXP && XLENGTH(held) > 0 &&
        static_cast<size_t>(XLENGTH(held)) < kRecordBytes) {
      have = static_cast<size_t>(XLENGTH(held));
      std::memcpy(bytes, RAW(held), have);
    }
  }

  bool take(Record& rec) {
    if (have < kRecordBytes) have += emc_pool::read_ready(rc, bytes + have, kRecordBytes - have);
    if (have < kRecordBytes) return false;
    int v[3];
    for (int k = 0; k < 3; ++k) {
      const unsigned char* b = bytes + 4 * k;
      v[k] = static_cast<int>(static_cast<unsigned int>(b[0]) |
                              (static_cast<unsigned int>(b[1]) << 8) |
                              (static_cast<unsigned int>(b[2]) << 16) |
                              (static_cast<unsigned int>(b[3]) << 24));
    }
    rec.w = v[0];
    rec.generation = v[1];
    rec.request = v[2];
    rec.status = bytes[12];
    have = 0;
    return true;
  }

  // The partial record back to the channel's buffer.
  void park() {
    if (TYPEOF(buf) != ENVSXP) return;
    RawVector held(have);
    if (have > 0) std::memcpy(RAW(held), bytes, have);
    set_var(buf, "bytes", held);
    have = 0;
  }
};

const int kStatusStart = 1;

enum Await { kRecord, kDead, kTimeout };

// .emc_wpool_await_record(): the next record of this pool generation that is
// not a start notice, which is noted in `picked`. A pending worker that
// exited gives kDead (with its index in rec.w); the limit passing with no
// pending worker computing gives kTimeout.
Await await_record(Records& records, Workers& ws, const std::vector<bool>& pending,
                   std::vector<bool>& picked, int generation, int request, double limit,
                   Record& rec, int& stale) {
  Steady::time_point started = Steady::now(), idle_since = started;
  long idle = 0;
  double wait = 0.0;
  bool waiting = false;
  for (;;) {
    if (records.take(rec)) {
      waiting = false;
      wait = 0.0;
      if (rec.generation != generation) {
        ++stale;
        continue;
      }
      if (rec.status == kStatusStart) {
        if (rec.request != request) ++stale;
        else if (rec.w >= 1 && rec.w <= ws.n) picked[rec.w - 1] = true;
        continue;
      }
      return kRecord;
    }
    ++idle;
    if ((idle % 64) == 0) {
      std::vector<int> open;
      for (int w = 0; w < ws.n; ++w) if (pending[w]) open.push_back(w);
      for (int w : open) {
        if (!ws.alive(w)) {
          rec.w = w + 1;
          return kDead;
        }
      }
      if (seconds_since(started) > limit) {
        if (!ws.busy(open)) return kTimeout;
        started = Steady::now();
      }
    }
    if (!waiting) {
      waiting = true;
      idle_since = Steady::now();
    }
    emc_pool::idle_wait(idle_since, wait);
    if (wait > 0 && (idle % 256) == 0) Rcpp::checkUserInterrupt();
  }
}

}  // namespace

namespace emc_pool {

SEXP lpt_partition(const double* cost, int n, int n_workers) {
  if (n_workers <= 1) {
    IntegerVector all(n);
    for (int s = 0; s < n; ++s) all[s] = s + 1;
    return List::create(all);
  }
  std::vector<int> order(n);
  for (int s = 0; s < n; ++s) order[s] = s;
  std::stable_sort(order.begin(), order.end(),
                   [cost](int a, int b) { return cost[a] > cost[b]; });
  std::vector<std::vector<int> > held(n_workers);
  std::vector<double> load(n_workers, 0.0);
  for (int s : order) {
    int w = 0;
    for (int k = 1; k < n_workers; ++k) if (load[k] < load[w]) w = k;
    held[w].push_back(s + 1);
    load[w] += cost[s];
  }
  List out(n_workers);
  for (int w = 0; w < n_workers; ++w) out[w] = IntegerVector(held[w].begin(), held[w].end());
  IntegerVector ord(n);
  for (int k = 0; k < n; ++k) ord[k] = order[k] + 1;
  out.attr("order") = ord;
  return out;
}

double cores_now(SEXP core_ctl, double base) {
  if (core_ctl == R_NilValue) return base;
  const int n_chains = Rf_asInteger(lget(core_ctl, "n_chains"));
  const int n_slots = Rf_asInteger(lget(core_ctl, "n_slots"));
  const int total = Rf_asInteger(lget(core_ctl, "total"));
  int finished = 0;
#ifndef _WIN32
  SEXP dir = lget(core_ctl, "dir");
  if (TYPEOF(dir) == STRSXP && XLENGTH(dir) == 1) {
    DIR* d = opendir(CHAR(STRING_ELT(dir, 0)));
    if (d != NULL) {
      while (struct dirent* e = readdir(d)) if (e->d_name[0] != '.') ++finished;
      closedir(d);
    }
  }
#else
  static Function now = package_function(".emc_cores_now");
  return Rf_asReal(now(core_ctl, base));
#endif
  const int active = std::min(n_slots, n_chains - finished);
  if (active <= 0) return total;
  return std::max(base, static_cast<double>(total / active));
}

void lost(RObject& pool, const char* why) {
  degraded(why);
  pool = lset(pool, "alive", Rf_ScalarLogical(FALSE));
}

List requests(RObject& pool, const List& msgs) {
  Workers ws(pool);
  const int n = std::min(ws.n, static_cast<int>(msgs.size()));
  List replies(msgs.size());
  std::vector<bool> sent(n, false);
  std::vector<unsigned char> payload;
  for (int w = 0; w < n; ++w) {
    if (Rf_isNull(msgs[w])) continue;
    serialize(msgs[w], payload);
    if (!ws.send(w, payload, deadline()).empty()) {
      lost(pool, transport_error("send to a worker failed").c_str());
      break;
    }
    sent[w] = true;
  }
  const double limit = deadline();
  for (int w = 0; w < n; ++w) if (sent[w]) replies[w] = ws.reply(w, limit);
  return replies;
}

bool subject_ll_queue(RObject& pool, SEXP rows, SEXP ctx, SEXP part, const char* source, int k,
                      double* ll) {
  Claim claim;
  if (TYPEOF(rows) != VECSXP || !claim_attach(pool, claim)) return false;
  const int n = static_cast<int>(XLENGTH(rows));
  Workers ws(pool);
  std::vector<int> order;
  if (TYPEOF(part) != VECSXP || XLENGTH(part) != ws.n || !queue_order(part, n, order))
    return false;
  SEXP data = lget(ctx, "data");
  RObject pooled;
  if (TYPEOF(data) == VECSXP && XLENGTH(data) > 0 && Rf_inherits(VECTOR_ELT(data, 0), "data.frame"))
    pooled = Rf_getAttrib(VECTOR_ELT(data, 0), Rf_install("pooled_values"));

  std::vector<unsigned char> bytes;
  {
    List shared = List::create(_["rows"] = rows,
                               _["order"] = IntegerVector(order.begin(), order.end()),
                               _["pooled"] = pooled);
    serialize(shared, bytes);
  }
  SharedFile file;
  const int request = next_request();
  std::vector<bool> sent(ws.n, false);
  bool alive = true;
  const std::string unpublished = shared_write(lget(pool, "dir"), bytes, file.path);
  if (!unpublished.empty()) {
    transport_fail(unpublished);
    degraded(unpublished);
    alive = false;
  } else {
    emc_claim_reset(claim.ptr, request);
    List msg = List::create(_["kind"] = "subject_ll_claim", _["shared_file"] = file.path,
                            _["claim_path"] = claim.path, _["key"] = request,
                            _["n"] = static_cast<int>(order.size()), _["notify"] = false);
    serialize(msg, bytes);
    const double limit = deadline();
    for (int w = 0; w < ws.n; ++w) {
      if (!ws.send(w, bytes, limit).empty()) {
        degraded(transport_error("send to a worker failed"));
        alive = false;
        break;
      }
      sent[w] = true;
    }
  }

  std::fill(ll, ll + static_cast<R_xlen_t>(n) * k, 0.0);
  std::vector<bool> stored(n, false);
  double worker_cpu = 0.0;
  List failures;
  auto keep_failures = [&failures](SEXP more) {
    if (TYPEOF(more) != VECSXP) return;
    for (R_xlen_t e = 0; e < XLENGTH(more); ++e) failures.push_back(VECTOR_ELT(more, e));
  };
  // A subject's values are numbers, one per candidate, none of them +Inf.
  auto values_ok = [k](SEXP x) {
    if (TYPEOF(x) != REALSXP || XLENGTH(x) != k) return false;
    for (int c = 0; c < k; ++c) if (ISNAN(REAL(x)[c]) || REAL(x)[c] == R_PosInf) return false;
    return true;
  };
  const double limit = deadline();
  std::vector<int> subs;
  for (int w = 0; w < ws.n; ++w) {
    if (!sent[w]) continue;
    RObject res = ws.reply(w, limit);
    SEXP vals = lget(res, "ll");
    bool ok = TYPEOF(res) == VECSXP && lget(res, "failed") == R_NilValue &&
        reply_subjects(lget(res, "subs"), n, stored, subs) && TYPEOF(vals) == VECSXP &&
        XLENGTH(vals) == static_cast<R_xlen_t>(subs.size());
    for (size_t e = 0; ok && e < subs.size(); ++e) ok = values_ok(VECTOR_ELT(vals, e));
    if (!ok) {
      emc_claim_reset(claim.ptr, 0);
      degraded("subject likelihood worker gave no valid reply");
      alive = false;
      continue;
    }
    for (size_t e = 0; e < subs.size(); ++e) {
      const double* x = REAL(VECTOR_ELT(vals, e));
      for (int c = 0; c < k; ++c) ll[subs[e] + static_cast<R_xlen_t>(c) * n] = x[c];
      stored[subs[e]] = true;
    }
    keep_failures(lget(res, "failures"));
    SEXP cpu = lget(res, "cpu");
    if ((TYPEOF(cpu) == REALSXP || TYPEOF(cpu) == INTSXP) && Rf_length(cpu) == 1)
      worker_cpu += Rf_asReal(cpu);
  }
  std::vector<int> missing;
  for (int s = 0; s < n; ++s) if (!stored[s]) missing.push_back(s + 1);
  if (!missing.empty()) {
    static Function compute = package_function(".emc_wpool_compute");
    IntegerVector miss(missing.begin(), missing.end());
    List rows_m(missing.size());
    for (size_t e = 0; e < missing.size(); ++e) rows_m[e] = VECTOR_ELT(rows, missing[e] - 1);
    List res = compute(List::create(_["kind"] = "subject_ll", _["subs"] = miss,
                                    _["rows"] = rows_m, _["pooled"] = pooled), ctx);
    const List vals = res["ll"];
    if (vals.size() != static_cast<R_xlen_t>(missing.size()))
      Rcpp::stop("subject likelihoods: the round's values differ in shape");
    for (size_t e = 0; e < missing.size(); ++e) {
      const NumericVector x = as<NumericVector>(vals[e]);
      if (x.size() != k) Rcpp::stop("subject likelihoods: the round's values differ in shape");
      for (int c = 0; c < k; ++c) ll[missing[e] - 1 + static_cast<R_xlen_t>(c) * n] = x[c];
    }
    keep_failures(res["failures"]);
  }
  if (failures.size() > 0) {
    static Function note = package_function(".emc_wpool_note_failures");
    note(failures, std::string(source));
  }
  state_add_real("subject_ll_cpu", worker_cpu);
  if (!alive) pool = lset(pool, "alive", Rf_ScalarLogical(FALSE));
  return true;
}

bool particle_queue(RObject& pool, const ParticleRequest& q, ParticleReply& out) {
  SEXP done = lget(pool, "done");
  SEXP profile = Rf_GetOption1(Rf_install("emc2.sampler_profile"));
  Claim claim;
  if (done == R_NilValue || flag(profile) || TYPEOF(q.pm_settings) != VECSXP ||
      !claim_attach(pool, claim))
    return false;
  const int n = static_cast<int>(XLENGTH(q.pm_settings)), p = q.n_pars;
  Workers ws(pool);
  std::vector<int> order;
  if (!queue_order(q.part, n, order)) return false;
  SEXP zero = emc_iter::package_value(".EMC_REJECT_ZERO");

  std::vector<unsigned char> bytes;
  {
    List inputs = List::create(
        _["alpha"] = q.alpha, _["population_mu"] = q.population_mu, _["pm"] = q.pm_settings,
        _["prev_ll"] = q.prev_ll, _["seeds"] = q.seeds, _["extra"] = q.extra,
        _["pooled"] = q.pooled, _["rows_only"] = q.rows_only, _["skip"] = q.skip,
        _["order"] = IntegerVector(order.begin(), order.end()));
    List shared = List::create(
        _["group_var"] = lget(q.group_chol, "ref"), _["idx_list"] = lget(q.group_chol, "idx_list"),
        _["marginal_idx"] = lget(q.group_chol, "marginal_idx"),
        _["include_full"] = lget(q.group_chol, "include_full"), _["inputs"] = inputs);
    serialize(shared, bytes);
  }
  SharedFile file;
  bool alive = true;
  const std::string unpublished = shared_write(lget(pool, "dir"), bytes, file.path);
  if (!unpublished.empty()) {
    transport_fail(unpublished);
    degraded(unpublished);
    alive = false;
  }
  SEXP gen = lget(pool, "generation");
  int generation = 0;
  if ((TYPEOF(gen) == INTSXP || TYPEOF(gen) == REALSXP) && Rf_length(gen) >= 1 &&
      Rf_asInteger(gen) != NA_INTEGER)
    generation = Rf_asInteger(gen);
  const int request = next_request();
  std::vector<bool> sent(ws.n, false);
  if (alive) {
    emc_claim_reset(claim.ptr, request);
    const double limit = deadline();
    for (int w = 0; w < ws.n; ++w) {
      List msg = List::create(_["kind"] = "claim", _["shared_file"] = file.path,
                              _["claim_path"] = claim.path, _["key"] = request,
                              _["n"] = static_cast<int>(order.size()), _["notify"] = true,
                              _["w"] = w + 1, _["generation"] = generation,
                              _["request"] = request);
      serialize(msg, bytes);
      if (ws.send(w, bytes, limit).empty()) {
        sent[w] = true;
      } else {
        degraded(transport_error("send to a worker failed"));
        alive = false;
      }
    }
  }
  const Steady::time_point receive_started = Steady::now();

  NumericMatrix props(p + 1, n);
  NumericVector times(n);
  List pm(Rf_shallow_duplicate(q.pm_settings)), seeds(Rf_shallow_duplicate(q.seeds));
  const bool with_extra = q.extra != R_NilValue;
  List extra_ll(with_extra ? n : 0);
  RObject rejects = Rf_duplicate(zero);
  RObject ens_counts = Rf_duplicate(zero);
  RObject ens_msgs = Rf_allocVector(VECSXP, 0);
  double cpu_total = 0.0;
  std::vector<bool> stored(n, false);
  std::vector<int> subs;

  // A reply's results into the round's, by the subjects it names. FALSE,
  // with nothing stored, when the reply does not fit them.
  auto store = [&](SEXP res, bool from_worker) {
    if (TYPEOF(res) != VECSXP || !reply_subjects(lget(res, "subs"), n, stored, subs)) return false;
    const R_xlen_t m = static_cast<R_xlen_t>(subs.size());
    SEXP r_pm = lget(res, "pm"), r_seeds = lget(res, "seeds");
    // Whole numbers arrive as integers.
    Shield<SEXP> r_props(TYPEOF(lget(res, "props")) == INTSXP
                             ? Rf_coerceVector(lget(res, "props"), REALSXP) : lget(res, "props"));
    Shield<SEXP> r_times(TYPEOF(lget(res, "times")) == INTSXP
                             ? Rf_coerceVector(lget(res, "times"), REALSXP) : lget(res, "times"));
    if (TYPEOF(r_props) != REALSXP || !Rf_isMatrix(r_props) || Rf_ncols(r_props) != m ||
        (m > 0 && Rf_nrows(r_props) != p + 1) || TYPEOF(r_pm) != VECSXP || XLENGTH(r_pm) != m ||
        TYPEOF(r_seeds) != VECSXP || XLENGTH(r_seeds) != m || TYPEOF(r_times) != REALSXP ||
        XLENGTH(r_times) != m)
      return false;
    SEXP r_extra = lget(res, "extra_ll");
    const bool has_extra = with_extra && TYPEOF(r_extra) == VECSXP && XLENGTH(r_extra) == m;
    for (R_xlen_t e = 0; e < m; ++e) {
      const int s = subs[e];
      std::copy(REAL(r_props) + e * (p + 1), REAL(r_props) + (e + 1) * (p + 1),
                props.begin() + static_cast<R_xlen_t>(s) * (p + 1));
      times[s] = REAL(r_times)[e];
      pm[s] = VECTOR_ELT(r_pm, e);
      seeds[s] = VECTOR_ELT(r_seeds, e);
      if (has_extra) extra_ll[s] = VECTOR_ELT(r_extra, e);
      stored[s] = true;
    }
    SEXP r_rej = lget(res, "rejects");
    if (r_rej != R_NilValue) rejects = add_counts(rejects, r_rej);
    SEXP r_ens = lget(res, "ensemble_rejects");
    if (r_ens != R_NilValue) {
      ens_counts = add_counts(ens_counts, lget(r_ens, "counts"));
      SEXP msgs = lget(r_ens, "msgs");
      SEXP names = Rf_getAttrib(msgs, R_NamesSymbol);
      for (R_xlen_t e = 0; TYPEOF(msgs) == VECSXP && e < Rf_xlength(names); ++e) {
        const char* cls = CHAR(STRING_ELT(names, e));
        if (lget(ens_msgs, cls) == R_NilValue && VECTOR_ELT(msgs, e) != R_NilValue)
          ens_msgs = lset(ens_msgs, cls, VECTOR_ELT(msgs, e));
      }
    }
    SEXP cpu = lget(res, "cpu");
    if (from_worker && (TYPEOF(cpu) == REALSXP || TYPEOF(cpu) == INTSXP) && Rf_length(cpu) == 1)
      cpu_total += Rf_asReal(cpu);
    return true;
  };

  std::vector<bool> pending(sent), picked(ws.n, false);
  int stale = 0;
  Records records(done);
  auto any = [](const std::vector<bool>& x) {
    return std::find(x.begin(), x.end(), true) != x.end();
  };
  if (any(pending)) {
    const double limit = deadline();
    while (any(pending)) {
      Record rec = {0, 0, 0, 0};
      const Await code = await_record(records, ws, pending, picked, generation, request, limit,
                                      rec, stale);
      if (code == kRecord && rec.request != request) {
        ++stale;
        continue;
      }
      const bool placed = code == kRecord && rec.w >= 1 && rec.w <= ws.n && pending[rec.w - 1];
      if (!placed) {
        emc_claim_reset(claim.ptr, 0);
        records.park();
        if (code == kTimeout) {
          IntegerVector stuck, idle;
          for (int w = 0; w < ws.n; ++w) if (pending[w]) {
            stuck.push_back(w + 1);
            if (!picked[w]) idle.push_back(w + 1);
          }
          static Function timed_out = package_function(".emc_wpool_iter_timeout");
          timed_out(pool, stuck, idle, limit);
        } else if (code == kDead) {
          degraded("worker gave no reply");
        } else {
          degraded("lost alignment on the worker completion channel");
        }
        alive = false;
        static Function retire = package_function(".emc_wpool_retire");
        retire(pool);
        break;
      }
      const int w = rec.w - 1;
      RObject res = ws.reply(w, limit);
      pending[w] = false;
      if (res.isNULL()) {
        emc_claim_reset(claim.ptr, 0);
        degraded(transport_error("worker gave no reply"));
        alive = false;
      } else if (lget(res, "failed") == R_NilValue) {
        store(res, true);
      }
    }
  }
  records.park();
  if (stale > 0) state_add_int("stale_records", stale);

  // Subjects no worker answered for run here.
  for (int s = 0; s < n; ++s) {
    if (stored[s]) continue;
    static Function compute = package_function(".emc_wpool_iter_missing");
    RObject res = compute(s + 1, q.ctx, q.alpha, q.population_mu, pm, q.prev_ll, seeds, q.extra,
                          q.pooled, q.skip, q.group_chol, q.rows_only);
    RObject named = lset(res, "subs", Rf_ScalarInteger(s + 1));
    if (!store(named, false))
      Rcpp::stop("particle step: a subject's result differs in shape");
  }
  record_receive(seconds_since(receive_started));

  out.props = props;
  out.pm_settings = pm;
  out.seeds = seeds;
  out.times = times;
  out.extra_ll = with_extra ? static_cast<SEXP>(extra_ll) : R_NilValue;
  out.rejects = rejects;
  out.ensemble_rejects = List::create(_["counts"] = ens_counts, _["msgs"] = ens_msgs);
  out.worker_cpu = cpu_total;
  if (!alive) pool = lset(pool, "alive", Rf_ScalarLogical(FALSE));
  return true;
}

}  // namespace emc_pool

// .emc_wpool_subject_ll() over the shared subject queue: the likelihoods (a
// vector for one candidate row) and the pool as the round left it. NULL when
// the pool does not serve the queue.
// [[Rcpp::export(rng = false)]]
SEXP emc_wpool_subject_ll_queue_cpp(SEXP pool, List rows, SEXP ctx, SEXP part,
                                    std::string source, SEXP ns) {
  emc_iter::set_package_namespace(ns);
  const int n = rows.size();
  if (n == 0 || !Rf_isMatrix(rows[0])) return R_NilValue;
  const int k = Rf_nrows(rows[0]);
  NumericMatrix ll(n, k);
  RObject after = pool;
  if (!emc_pool::subject_ll_queue(after, rows, ctx, part, source.c_str(), k, ll.begin()))
    return R_NilValue;
  if (k == 1) ll.attr("dim") = R_NilValue;
  return List::create(_["ll"] = ll, _["pool"] = after);
}

// .emc_wpool_iter() over the shared subject queue, unprofiled. NULL when the
// pool does not serve the queue.
// [[Rcpp::export(rng = false)]]
SEXP emc_wpool_particle_queue_cpp(SEXP pool, SEXP ctx, SEXP part, SEXP alpha,
                                  SEXP population_mu, SEXP group_chol, SEXP pm_settings,
                                  SEXP prev_ll, SEXP seeds, SEXP extra, SEXP pooled, SEXP skip,
                                  int n_pars, bool rows_only, SEXP ns) {
  emc_iter::set_package_namespace(ns);
  emc_pool::ParticleRequest request;
  request.ctx = ctx;
  request.part = part;
  request.alpha = alpha;
  request.population_mu = population_mu;
  request.group_chol = group_chol;
  request.pm_settings = pm_settings;
  request.prev_ll = prev_ll;
  request.seeds = seeds;
  request.extra = extra;
  request.pooled = pooled;
  request.skip = skip;
  request.n_pars = n_pars;
  request.rows_only = rows_only;
  emc_pool::ParticleReply reply;
  RObject after = pool;
  if (TYPEOF(seeds) != VECSXP || TYPEOF(pm_settings) != VECSXP ||
      XLENGTH(seeds) != XLENGTH(pm_settings) ||
      !emc_pool::particle_queue(after, request, reply))
    return R_NilValue;
  return List::create(_["props"] = reply.props, _["pm_settings"] = reply.pm_settings,
                      _["seeds"] = reply.seeds, _["times"] = reply.times,
                      _["extra_ll"] = reply.extra_ll, _["alive"] = lget(after, "alive"),
                      _["rejects"] = reply.rejects,
                      _["ensemble_rejects"] = reply.ensemble_rejects,
                      _["worker_cpu"] = reply.worker_cpu, _["profile"] = R_NilValue);
}
