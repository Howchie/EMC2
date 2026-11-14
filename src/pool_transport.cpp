// The worker pool's byte transport (R/chain_pool.R): messages as serialized
// payloads in length-prefixed segments over the pool's connections. A segment
// is an 8-byte little-endian double, its length (negated when more segments
// follow), then that many bytes. Reading polls a non-blocking connection,
// backing off, until the bytes, the deadline or the other end's exit.
#include <Rcpp.h>
// The connection structure is a C declaration with members named `class`
// and `private`.
#define class connection_class
#define private connection_private
extern "C" {
#include <R_ext/Connections.h>
}
#undef class
#undef private
#include "pool_transport.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

using namespace Rcpp;

namespace {

typedef std::chrono::steady_clock Clock;

// How long a read polls without sleeping before it backs off.
const double kSpinSeconds = 1e-3;

double seconds_since(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

std::string format_seconds(const char* fmt, const char* what, double s) {
  char buf[160];
  if (what != NULL) std::snprintf(buf, sizeof buf, fmt, what, s);
  else std::snprintf(buf, sizeof buf, fmt, s);
  return std::string(buf);
}

bool is_true(SEXP x) {
  return TYPEOF(x) == LGLSXP && Rf_length(x) == 1 && LOGICAL(x)[0] == TRUE;
}

const double kReadChunk = 65536, kWriteChunk = 4096;
const double kMaxFrame = 2147483648.0, kMaxSegment = 1073741824.0;

// n bytes from the connection into `out`. An empty string on success, else
// what went wrong.
std::string get_bytes(Rconnection con, double n, const emc_pool::Watch& watch,
                      double read_chunk, std::vector<unsigned char>& out) {
  out.clear();
  if (n <= 0) return "";
  const size_t want_total = static_cast<size_t>(n);
  out.resize(want_total);
  const bool timed = watch.timed;
  const double limit = watch.limit;
  size_t got = 0;
  Clock::time_point started = Clock::now();
  double wait = 0.0;
  long idle = 0;
  bool waiting = false;
  Clock::time_point idle_since = started;
  while (got < want_total) {
    const size_t want = std::min(want_total - got, static_cast<size_t>(read_chunk));
    size_t k = 0;
    bool failed = false;
    if (!con->isopen || !con->canread) {
      if (!timed || !con->isopen) return "cannot read from this connection";
      failed = true;
    } else {
      k = R_ReadConnection(con, out.data() + got, want);
      // A read that reports more than was asked for failed (nothing ready on
      // a non-blocking pipe).
      if (k > want) {
        if (!timed) return "the read from the pipe failed";
        failed = true;
      }
    }
    if (failed) k = 0;
    if (k == 0) {
      if (!timed)
        return got == 0 ? "the other end of the pipe closed"
                        : "the pipe closed part-way through a message";
      ++idle;
      if ((idle % 64) == 0 && watch.alive) {
        if (!watch.alive())
          return got == 0 ? "the worker exited without replying"
                          : "the worker exited part-way through a message";
      }
      if (seconds_since(started) > limit) {
        if (watch.busy && watch.busy()) {
          started = Clock::now();
          continue;
        }
        return format_seconds("no %s within %.0f s", got == 0 ? "reply" : "rest of the reply",
                              limit);
      }
      // A reply that is about to arrive is taken at once: the wait polls
      // without sleeping at first, then backs off.
      if (!waiting) {
        waiting = true;
        idle_since = Clock::now();
      }
      emc_pool::idle_wait(idle_since, wait);
      if (wait > 0 && (idle % 256) == 0) Rcpp::checkUserInterrupt();
      continue;
    }
    waiting = false;
    wait = 0.0;
    got += k;
  }
  return "";
}

// The bytes to the connection in chunks; an empty string on success.
std::string put_bytes(Rconnection con, const unsigned char* bytes, size_t n, bool timed,
                      double limit, double write_chunk) {
  const Clock::time_point started = Clock::now();
  size_t from = 0;
  while (from < n) {
    if (timed && from > 0 && seconds_since(started) > limit)
      return format_seconds("the request was still being written after %.0f s", NULL, limit);
    const size_t len = std::min(n - from, static_cast<size_t>(write_chunk));
    if (!con->isopen || !con->canwrite) return "cannot write to this connection";
    const size_t wrote = R_WriteConnection(con, const_cast<unsigned char*>(bytes) + from, len);
    if (wrote != len) return "problem writing to connection";
    from += len;
  }
  return "";
}

SEXP failure(const std::string& msg) { return Rf_mkString(msg.c_str()); }

std::string recv_frames(Rconnection rc, const emc_pool::Watch& watch, double read_chunk,
                        double max_frame, std::vector<unsigned char>& payload) {
  std::vector<unsigned char> part, header;
  payload.clear();
  bool first = true;
  for (;;) {
    std::string err = get_bytes(rc, 8, watch, read_chunk, header);
    if (!err.empty()) return err;
    double n;
    std::memcpy(&n, header.data(), 8);
#if defined(WORDS_BIGENDIAN)
    {
      unsigned char* b = reinterpret_cast<unsigned char*>(&n);
      for (int i = 0; i < 4; ++i) std::swap(b[i], b[7 - i]);
    }
#endif
    const bool more = R_FINITE(n) && n < 0;
    if (more) n = -n;
    if (!R_FINITE(n) || n < 0 || n > max_frame) return "implausible message length on the pipe";
    err = get_bytes(rc, n, watch, read_chunk, part);
    if (!err.empty()) return err;
    if (first && !more) {
      payload.swap(part);
      break;
    }
    payload.insert(payload.end(), part.begin(), part.end());
    first = false;
    if (!more) break;
  }
  return "";
}

std::string send_frames(Rconnection wc, const unsigned char* bytes, size_t total, bool timed,
                        double limit, double max_segment, double write_chunk) {
  size_t from = 0;
  for (;;) {
    const size_t len = std::min(total - from, static_cast<size_t>(max_segment));
    const bool last = from + len >= total;
    double tag = last ? static_cast<double>(len) : -static_cast<double>(len);
    unsigned char header[8];
    std::memcpy(header, &tag, 8);
#if defined(WORDS_BIGENDIAN)
    for (int i = 0; i < 4; ++i) std::swap(header[i], header[7 - i]);
#endif
    std::string err = put_bytes(wc, header, 8, timed, limit, write_chunk);
    if (!err.empty()) return err;
    err = put_bytes(wc, bytes + from, len, timed, limit, write_chunk);
    if (!err.empty()) return err;
    if (last) break;
    from += len;
  }
  return "";
}

struct ConnectionLookup {
  SEXP object;
  Rconnection con;
};
void get_connection(void* data) {
  ConnectionLookup* look = static_cast<ConnectionLookup*>(data);
  look->con = R_GetConnection(look->object);
}

}  // namespace

namespace emc_pool {

void* connection(SEXP con) {
  if (!Rf_inherits(con, "connection")) return NULL;
  // The lookup raises an R error for a closed connection.
  ConnectionLookup look = {con, NULL};
  if (!R_ToplevelExec(get_connection, &look)) return NULL;
  return look.con;
}

std::string recv_payload(void* con, const Watch& watch, std::vector<unsigned char>& payload) {
  if (con == NULL) return "cannot read from this connection";
  return recv_frames(static_cast<Rconnection>(con), watch, kReadChunk, kMaxFrame, payload);
}

std::string send_payload(void* con, const unsigned char* bytes, size_t n, bool timed,
                         double limit) {
  if (con == NULL) return "cannot write to this connection";
  Rconnection wc = static_cast<Rconnection>(con);
  std::string err = send_frames(wc, bytes, n, timed, limit, kMaxSegment, kWriteChunk);
  if (err.empty() && wc->fflush != NULL) wc->fflush(wc);
  return err;
}

size_t read_ready(void* con, unsigned char* out, size_t n) {
  Rconnection rc = static_cast<Rconnection>(con);
  if (rc == NULL || !rc->isopen || !rc->canread) return 0;
  const size_t k = R_ReadConnection(rc, out, n);
  return k > n ? 0 : k;
}

void idle_wait(const Clock::time_point& idle_since, double& wait) {
  if (seconds_since(idle_since) <= kSpinSeconds) return;
  if (wait > 0) std::this_thread::sleep_for(std::chrono::duration<double>(wait));
  wait = std::min(std::max(wait * 2, 5e-5), 2e-3);
}

}  // namespace emc_pool

// One message's payload from the connection, its segments joined; a string
// saying what went wrong otherwise.
// [[Rcpp::export(rng = false)]]
SEXP emc_wpool_recv_cpp(SEXP con, SEXP deadline, SEXP alive, SEXP busy, double read_chunk,
                        double max_frame) {
  Rconnection rc = R_GetConnection(con);
  emc_pool::Watch watch;
  watch.timed = deadline != R_NilValue;
  watch.limit = watch.timed ? Rf_asReal(deadline) : 0.0;
  if (alive != R_NilValue) watch.alive = [alive]() { return is_true(Function(alive)()); };
  if (TYPEOF(busy) == CLOSXP) watch.busy = [busy]() { return is_true(Function(busy)()); };
  std::vector<unsigned char> payload;
  const std::string err = recv_frames(rc, watch, read_chunk, max_frame, payload);
  if (!err.empty()) return failure(err);
  RawVector out(payload.size());
  if (!payload.empty()) std::memcpy(RAW(out), payload.data(), payload.size());
  return out;
}

// A message's payload to the connection in segments; TRUE, or a string
// saying what went wrong.
// [[Rcpp::export(rng = false)]]
SEXP emc_wpool_send_cpp(SEXP con, RawVector payload, SEXP deadline, double max_segment,
                        double write_chunk) {
  Rconnection wc = R_GetConnection(con);
  const bool timed = deadline != R_NilValue;
  const std::string err = send_frames(wc, RAW(payload), static_cast<size_t>(XLENGTH(payload)),
                                      timed, timed ? Rf_asReal(deadline) : 0.0, max_segment,
                                      write_chunk);
  if (!err.empty()) return failure(err);
  return Rf_ScalarLogical(TRUE);
}
