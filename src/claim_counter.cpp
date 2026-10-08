// A shared work counter for the worker pool's subject queue.
//
// Workers claim subjects themselves, costliest first, instead of the chain
// process handing each one out: the chain then handles one message per worker
// per iteration rather than one per subject. The counter is a 64-bit word in a
// file mapped shared by every worker: the high half is the round's key and the
// low half the next index. A claim for any other key fails, so a worker still
// running from an earlier round cannot take work from the current one.

#include <Rcpp.h>
#include <cstdint>
#include <string>

#ifndef _WIN32
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace {

struct ClaimMap {
  int fd;
  std::uint64_t* word;
};

void claim_close(ClaimMap* m) {
  if (m == nullptr) return;
#ifndef _WIN32
  if (m->word != nullptr) munmap(static_cast<void*>(m->word), sizeof(std::uint64_t));
  if (m->fd >= 0) close(m->fd);
#endif
  m->word = nullptr;
  m->fd = -1;
}

void claim_finalize(SEXP p) {
  ClaimMap* m = static_cast<ClaimMap*>(R_ExternalPtrAddr(p));
  if (m == nullptr) return;
  claim_close(m);
  delete m;
  R_ClearExternalPtr(p);
}

ClaimMap* claim_get(SEXP p) {
  if (TYPEOF(p) != EXTPTRSXP) Rcpp::stop("not a claim counter");
  ClaimMap* m = static_cast<ClaimMap*>(R_ExternalPtrAddr(p));
  if (m == nullptr || m->word == nullptr) Rcpp::stop("claim counter is closed");
  return m;
}

}  // namespace

//' Map the claim counter at `path`, creating the file if needed
//' @noRd
// [[Rcpp::export]]
SEXP emc_claim_open(std::string path) {
#ifdef _WIN32
  Rcpp::stop("claim counters need POSIX shared mappings");
#else
  int fd = open(path.c_str(), O_RDWR | O_CREAT, 0600);
  if (fd < 0) Rcpp::stop("cannot open claim counter file");
  struct stat st;
  if (fstat(fd, &st) != 0) { close(fd); Rcpp::stop("cannot stat claim counter file"); }
  if (st.st_size < static_cast<off_t>(sizeof(std::uint64_t)) &&
      ftruncate(fd, sizeof(std::uint64_t)) != 0) {
    close(fd);
    Rcpp::stop("cannot size claim counter file");
  }
  void* addr = mmap(nullptr, sizeof(std::uint64_t), PROT_READ | PROT_WRITE,
                    MAP_SHARED, fd, 0);
  if (addr == MAP_FAILED) { close(fd); Rcpp::stop("cannot map claim counter file"); }
  ClaimMap* m = new ClaimMap{fd, static_cast<std::uint64_t*>(addr)};
  SEXP p = PROTECT(R_MakeExternalPtr(m, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(p, claim_finalize, TRUE);
  UNPROTECT(1);
  return p;
#endif
}

//' Start round `key`: the next claim returns index 0
//' @noRd
// [[Rcpp::export]]
void emc_claim_reset(SEXP map, int key) {
  ClaimMap* m = claim_get(map);
  const std::uint64_t word = static_cast<std::uint64_t>(static_cast<std::uint32_t>(key)) << 32;
  __atomic_store_n(m->word, word, __ATOMIC_SEQ_CST);
}

//' Claim the next index of round `key` below `n`; -1 once the round is
//' exhausted or has been replaced by another
//' @noRd
// [[Rcpp::export]]
int emc_claim_next(SEXP map, int key, int n) {
  ClaimMap* m = claim_get(map);
  const std::uint32_t want = static_cast<std::uint32_t>(key);
  std::uint64_t cur = __atomic_load_n(m->word, __ATOMIC_SEQ_CST);
  for (;;) {
    if (static_cast<std::uint32_t>(cur >> 32) != want) return -1;
    const std::uint32_t idx = static_cast<std::uint32_t>(cur & 0xffffffffu);
    if (n <= 0 || idx >= static_cast<std::uint32_t>(n)) return -1;
    if (__atomic_compare_exchange_n(m->word, &cur, cur + 1, false,
                                    __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)) {
      return static_cast<int>(idx);
    }
  }
}
