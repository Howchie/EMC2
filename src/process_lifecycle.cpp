// A worker must not outlive the fit that owns it.
//
// The clean backend already cascades: the master dies, its end of the
// template's request pipe closes, the template sees end-of-stream and kills
// the workers it forked.  Measured, including from SIGKILL.
//
// The fork backend has no such intermediary.  Its workers are the master's own
// `mcparallel` children, blocked on a request pipe, and an orderly shutdown is
// the only thing that ends them.  `kill -9` and `kill` on the master therefore
// left every worker running -- reparented to init, blocked forever, each
// holding a whole likelihood context.  That backend is not exotic: it is what
// `devtools::load_all()` selects, because the namespace is then not an
// installed package, and what any model with custom trend kernels selects,
// because external pointers cannot be serialised to a spawned process.
//
// PR_SET_PDEATHSIG is the kernel primitive for exactly this: the child asks to
// be signalled when its parent dies, whatever kills the parent and whether or
// not the parent gets to run any cleanup.
//
// It is set in the child after the fork, so there is a window in which the
// parent can die before the request is registered -- and a signal that was
// never armed will never arrive.  The check that closes that window is
// `getppid()`, which R does not expose, so both halves happen here rather than
// leaving the caller to do the interesting one.

#include <Rcpp.h>

#ifndef _WIN32
#include <unistd.h>
#endif

#ifdef __linux__
#include <sys/prctl.h>
#include <csignal>
#endif

//' Ask the kernel to signal this process when its parent dies
//'
//' @return 1 if armed and the parent is still there, -1 if armed but the
//'   parent had already gone (the caller should exit rather than wait for a
//'   signal that will not come), 0 where the platform has no equivalent and
//'   the caller should carry on as before.
//' @noRd
// [[Rcpp::export]]
int emc_arm_parent_death() {
#ifdef __linux__
  if (prctl(PR_SET_PDEATHSIG, SIGKILL) != 0) return 0;
  // Reparenting is the observable consequence of the parent's death, whatever
  // adopts the orphan -- init, or a subreaper.
  return (getppid() <= 1) ? -1 : 1;
#else
  return 0;
#endif
}
