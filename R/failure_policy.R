# What the sampler does when something fails, and why.
#
# This file is the policy, written down. It exists because the mechanism it
# governs has already misled us once: a diverging PDE solve inside the RLF/FPE
# corner wedged a sampler, and because every error on the particle path is
# caught and answered by repeating the previous state, the fit looked like fast,
# well-mixing sampling from the outside. Nothing was wrong with the catching --
# it is what makes the sampler robust -- but with the policy left implicit there
# was no way to tell a rejected proposal from a broken worker.
#
# Two places catch everything:
#
#   safe_new_particle()   (R/sampling.R) repeats the subject's previous state
#   the group Gibbs step  (R/sampling.R) repeats the previous iteration
#
# Repeating the previous state is not a workaround. Under *numerical* rejection
# it is the correct move: it is what preserves the PMwG invariant, because a
# proposal that lands where the model has no density must be rejected, and
# rejection means keeping the current state. Any policy here has to leave that
# untouched.
#
# The problem is that the same catch also swallows a failed allocation, a dead
# worker, a missing file and a genuine bug in a model's R code, and answers all
# of them the same way. C1 added the classification and the counters, which
# changed nothing. This file decides what to do with them.
#
# --- the classes ------------------------------------------------------------
#
# `.emc_classify_failure()` sorts a condition into one of four classes. The
# classifier is a heuristic over error messages, and that matters for the policy
# below: a rule that aborts a fit on a heuristic will eventually abort a fit it
# should not have.
#
#   numerical       the model has no density at this proposal, or a solver could
#                   not converge there. Out-of-bounds parameters, a
#                   non-positive-definite covariance, a non-finite likelihood.
#                   Expected, and rejection is the right answer.
#
#   infrastructure  the computation was never attempted: an allocation failure,
#                   a closed connection, a worker that died, a missing shared
#                   broadcast. Repeating the previous state here is not a
#                   rejection, it is a silently lost update -- the chain does
#                   not advance and nothing says so.
#
#   programming     an R error that indicates a defect: an object not found, an
#                   unused argument, a subscript out of bounds, non-conformable
#                   arguments. Repeating the previous state produces a
#                   plausible-looking posterior from code that never ran.
#
#   unknown         everything else. Treated as suspicious, not as numerical:
#                   assuming the benign class is how the first mistake happened.
#
# --- the policy -------------------------------------------------------------
#
# `options(emc2.failure_policy = ...)`, one of:
#
#   "report"  (default)  Nothing aborts. Every class is counted as before, and
#                        the first non-numerical failure of each kind raises an
#                        immediate warning naming the class, where it happened
#                        and the message; a stage that saw any of them prints a
#                        summary when it finishes. The fit behaves exactly as it
#                        did before this file existed -- what changes is that
#                        the failures are no longer invisible.
#
#   "strict"             Numerical failures still reject and continue. An
#                        infrastructure, programming or unknown failure is
#                        re-raised where it was caught, which ends the chain.
#                        For a user who would rather find out than get a
#                        posterior built partly from repeated states.
#
#   "silent"             The behaviour before C1: count, say nothing. For a
#                        long production run whose failure modes are already
#                        understood, and for reproducing an older result.
#
# The default is "report" and not "strict" deliberately. The audit's own warning
# applies: "a wide fit that previously limped would now stop". Aborting on a
# message-matching heuristic would turn a misclassified numerical corner into a
# lost overnight fit, which is a worse failure than the one being fixed. Making
# the failures visible is the part that is safe to do by default; deciding that
# they are fatal is the user's call.
#
# --- what happens mid-block -------------------------------------------------
#
# Under "report" and "silent", nothing changes: the iteration completes with a
# repeated state, the block finishes, and the samples up to that point are
# preserved and saved as usual.
#
# Under "strict", the error propagates out of `run_stage()`. In a worker, it is
# caught by the pool and returned as a failed reply; the master then recomputes
# that share itself, hits the same error, and re-raises it -- so a strict abort
# happens in the master, where the traceback is useful, and not in a forked
# process whose output may be lost. `run_emc()` receives the error from
# `auto_mclapply` for that chain. The last completed block has already been
# written to the checkpoint file if one was given, so a strict abort loses the
# current block, not the fit.
#
# --- what the user sees -----------------------------------------------------
#
# First occurrence, per class and per source, immediately:
#
#   EMC2: 3 particle updates were rejected for an infrastructure failure, not a
#   numerical one ... the chain does not advance ... first message: <message>
#
# and at the end of a stage in which any non-numerical failure was counted, one
# summary line per class with its count and its first message.
#
# The warning is `immediate. = TRUE` and still best-effort, because of where it
# is raised: a particle failure is classified inside the worker, and a forked
# worker's stderr reaches the master while a *spawned* worker's goes to the
# template's log file, where nobody is reading. The counts travel back to the
# master with the reply, so the stage summary is the channel that always
# arrives; the warning is what makes a failure visible in the same second it
# happens, when it can be.
#
# One gap, stated rather than hidden: the stage summary is assembled in the
# master from counts the workers return with their replies. Where no worker
# pool can be started at all and there is more than one subject -- Windows, or a
# machine without `mkfifo` -- the particle step falls back to `mcmapply`, whose
# children have no channel back except the proposal matrix, and their counts are
# lost with the process. The immediate warning is still raised there, which is
# why it is `immediate. = TRUE` rather than deferred, but it is best-effort.
# Every route that uses the pool, including the single-worker route that
# computes in the master, reports in full.
#
# Numerical rejections are not announced. They are ordinary, they can be
# frequent in a hard corner of a model, and a warning per occurrence would train
# the user to ignore the channel that carries the other three.

.EMC_FAILURE_POLICIES <- c("silent", "report", "strict")

# Classes that are never allowed to abort, whatever the policy says. Rejecting a
# proposal the model has no density for is the sampler working.
.EMC_FAILURE_BENIGN <- "numerical"

.emc_failure_policy <- function() {
  opt <- getOption("emc2.failure_policy", "report")
  if (length(opt) != 1L || !is.character(opt) ||
      !(opt %in% .EMC_FAILURE_POLICIES)) {
    return("report")
  }
  opt
}

# "continue" repeats the previous state, which is what every caller did before
# this file existed. "abort" re-raises at the catch site.
.emc_failure_action <- function(cls) {
  if (identical(cls, .EMC_FAILURE_BENIGN)) return("continue")
  if (identical(.emc_failure_policy(), "strict")) "abort" else "continue"
}

# One warning per (source, class) per process. A wide fit can reject thousands
# of times; the count is what says how often, and the message is what says what.
.emc_failure_announce <- function(cls, source, cond) {
  if (identical(cls, .EMC_FAILURE_BENIGN)) return(invisible(FALSE))
  if (identical(.emc_failure_policy(), "silent")) return(invisible(FALSE))
  key <- paste(source, cls, sep = ".")
  seen <- .emc_profile_state$announced
  if (is.null(seen)) seen <- character(0)
  if (key %in% seen) return(invisible(FALSE))
  .emc_profile_state$announced <- c(seen, key)
  msg <- tryCatch(conditionMessage(cond), error = function(e) "")
  if (!length(msg) || !nzchar(msg[[1L]])) msg <- "(no message)"
  # `immediate. = TRUE`: this may be running in a forked worker, where a
  # deferred warning goes with the process.
  warning(sprintf(
    paste0("EMC2: a %s update failed for an %s reason, not a numerical one. ",
           "The previous state was repeated, so this chain did not advance ",
           "at that point and the fit will not say so on its own. ",
           "First message: %s\n",
           "  options(emc2.failure_policy = \"strict\") to stop instead."),
    source, cls, msg[[1L]]),
    call. = FALSE, immediate. = TRUE)
  invisible(TRUE)
}

# Re-raise at the catch site, with the class attached so a handler further out
# can tell why the fit stopped.
.emc_failure_abort <- function(cls, source, cond) {
  msg <- tryCatch(conditionMessage(cond), error = function(e) "")
  if (!length(msg) || !nzchar(msg[[1L]])) msg <- "(no message)"
  stop(errorCondition(
    sprintf(paste0("EMC2 stopped on a %s failure in the %s update ",
                   "(emc2.failure_policy = \"strict\"): %s"),
            cls, source, msg[[1L]]),
    class = c("emc_failure", paste0("emc_failure_", cls)),
    emc_class = cls, emc_source = source, emc_cause = cond))
}

# --- stage summary ----------------------------------------------------------

# Totals for one stage, as a data frame of the non-numerical classes that
# actually occurred. `counts` is a named list of per-source count vectors, the
# shape `.emc_reject_counts()` returns.
.emc_failure_summary <- function(counts) {
  rows <- list()
  for (source in names(counts)) {
    cur <- counts[[source]]
    if (is.null(cur)) next
    for (cls in names(cur)) {
      if (identical(cls, .EMC_FAILURE_BENIGN) || !isTRUE(cur[[cls]] > 0L)) next
      key <- paste(source, cls, sep = ".")
      msgs <- .emc_reject_messages()
      rows[[length(rows) + 1L]] <- data.frame(
        source = source, class = cls, n = as.integer(cur[[cls]]),
        message = if (is.null(msgs[[key]])) "" else msgs[[key]],
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(-out$n), , drop = FALSE]
}

# Emitted once per stage, in the master, where a warning cannot be lost with a
# forked process.  Silent when the only failures were numerical, which is the
# ordinary case and not news.
.emc_failure_report_stage <- function(counts, stage = NULL) {
  if (identical(.emc_failure_policy(), "silent")) return(invisible(NULL))
  summary <- .emc_failure_summary(counts)
  if (is.null(summary)) return(invisible(NULL))
  lines <- vapply(seq_len(nrow(summary)), function(i) {
    sprintf("  %d %s update%s classified as %s: %s", summary$n[i],
            summary$source[i], if (summary$n[i] > 1L) "s" else "",
            summary$class[i], summary$message[i])
  }, character(1))
  warning(sprintf(
    paste0("EMC2%s: %d update%s were answered by repeating the previous state ",
           "for a reason that was not numerical rejection. The chain did not ",
           "advance at those points.\n%s"),
    if (is.null(stage)) "" else paste0(" (", stage, " stage)"),
    sum(summary$n), if (sum(summary$n) > 1L) "s" else "",
    paste(lines, collapse = "\n")),
    call. = FALSE, immediate. = TRUE)
  invisible(summary)
}
