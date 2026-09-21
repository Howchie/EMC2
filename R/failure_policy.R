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
# The particle path catches numerical proposal failures and repeats the
# subject's previous state.  A group Gibbs draw is different: it is the state
# transition itself, so a failed draw aborts the chain with its location and
# numerical diagnostics rather than manufacturing a repeated iteration.
#
# Repeating a subject state is not a workaround. Under *particle numerical*
# rejection it is the correct move: a proposal that lands where the model has
# no density must be rejected, and rejection means keeping the current state.
#
# The particle catch can also swallow a failed allocation, a dead worker, a
# missing file and a genuine bug in a model's R code, and answer all of them
# the same way. C1 added the classification and counters; this file decides
# what to do with those particle-path failures.
#
# --- the classes ------------------------------------------------------------
#
# `.emc_classify_failure()` sorts a condition into one of four classes. The
# classifier is a heuristic over error messages, and that matters for the policy
# below: a rule that aborts a fit on a heuristic will eventually abort a fit it
# should not have.
#
#   numerical       the model has no density at this particle proposal, or a
#                   solver could not converge there. Out-of-bounds parameters,
#                   a non-positive-definite covariance, a non-finite
#                   likelihood. Expected on the particle path; a group Gibbs
#                   failure is fatal regardless of this classification.
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
#   "report"  (default)  Particle numerical failures continue as rejections.
#                        Every non-numerical particle failure is counted and
#                        announced; group Gibbs failures always abort with
#                        context. A stage that saw non-numerical particle
#                        failures prints a summary when it finishes.
#
#   "strict"             Same Gibbs behavior. On the particle path, numerical
#                        failures still reject and continue; infrastructure,
#                        programming and unknown failures are re-raised.
#
#   "silent"             The behaviour before C1: count, say nothing. For a
#                        long production run whose failure modes are already
#                        understood, and for reproducing an older result.
#
# The default is "report" for particle failures. Group Gibbs errors are not
# proposal rejections and are therefore fatal under every policy.
#
# --- what happens mid-block -------------------------------------------------
#
# Under "report" and "silent", particle numerical rejections still complete
# with a repeated subject state. A failed group Gibbs draw aborts immediately;
# it is never converted into a complete repeated iteration.
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

# A group Gibbs draw is a state transition, not a proposal.  Repeating the
# previous complete iteration after it fails therefore hides a broken chain
# rather than preserving a valid rejection.  Keep this separate from the
# particle rejection policy: numerical particle failures may continue, while a
# failed group update is always fatal and carries its stage/iteration context.
.emc_gibbs_abort <- function(cond, stage = NULL, iteration = NULL,
                             nuisance = FALSE) {
  msg <- tryCatch(conditionMessage(cond), error = function(e) "")
  if (!length(msg) || !nzchar(msg[[1L]])) msg <- "(no message)"
  where <- if (isTRUE(nuisance)) "nuisance group Gibbs" else "group Gibbs"
  location <- paste0(
    if (!is.null(stage)) paste0(" stage=", stage) else "",
    if (!is.null(iteration)) paste0(" iteration=", iteration) else ""
  )
  stop(errorCondition(
    sprintf("EMC2 aborted after a failed %s update%s: %s",
            where, location, msg[[1L]]),
    class = c("emc_gibbs_failure", "emc_failure", "emc_failure_numerical"),
    emc_class = "numerical", emc_source = "gibbs",
    emc_stage = stage, emc_iteration = iteration,
    emc_nuisance = isTRUE(nuisance), emc_cause = cond))
}

# A forked particle worker that dies BELOW R level -- a segfault, or the OOM
# killer -- never raises an R error, so `safe_new_particle`'s tryCatch cannot
# see it.  `parallel::mclapply` then simply DROPS that element: the result comes
# back short, with no indication of which subject is missing.  Assembling it
# positionally therefore slid the survivors into the wrong subjects' slots, and
# unlisting into a fixed-size array recycled a survivor into the gap -- one dead
# worker silently overwrote another subject's draw and the chain wedged with
# nothing reported.  `.emc_assemble_proposals()` matches results to subjects by
# tag instead; this is the channel that says a worker was lost.
.emc_worker_loss <- function(lost, stage = NULL, iteration = NULL) {
  location <- paste0(
    if (!is.null(stage)) paste0(" stage=", stage) else "",
    if (!is.null(iteration)) paste0(" iteration=", iteration) else ""
  )
  subs <- paste(lost, collapse = ", ")
  if (identical(.emc_failure_policy(), "strict")) {
    stop(errorCondition(
      sprintf(paste0("EMC2 lost %d particle worker(s)%s (subject(s) %s). The ",
                     "worker died below R level, so no R error was raised."),
              length(lost), location, subs),
      class = c("emc_worker_lost", "emc_failure"),
      emc_class = "crashed", emc_source = "particle",
      emc_stage = stage, emc_iteration = iteration, emc_subjects = lost))
  }
  if (identical(.emc_failure_policy(), "silent")) return(invisible(FALSE))
  # Not once-per-process like `.emc_failure_announce`: a lost worker is rare and
  # each loss costs a real update, so every one is worth a line.
  warning(sprintf(
    paste0("EMC2: %d particle worker(s) died%s (subject(s) %s). Those subjects ",
           "repeat their previous state for this iteration; every other ",
           "subject is unaffected.\n",
           "  options(emc2.failure_policy = \"strict\") to stop instead."),
    length(lost), location, subs), call. = FALSE, immediate. = TRUE)
  invisible(TRUE)
}

# How many consecutive TRIES a chain may go without moving at all before the
# fit is stopped -- the "[stage | try=N | ...]" blocks the console reports,
# step_size iterations each.  Not iterations: burn-in routinely rejects for
# tens of iterations in a row and then moves on (RDMGBM burn reached 32), and a
# limit of 25 iterations killed RDM, RDMGBM and LBA fits that were about to
# recover.  A healthy chain moves many times within one try -- 0/12 healthy
# fits tripped even a ONE-try limit -- so three whole tries without a single
# accepted move only happens to a chain that is genuinely stuck.
# Inf disables the check.
.emc_stall_try_limit <- function() {
  lim <- getOption("emc2.stall_try_limit", 3L)
  if (length(lim) != 1L || !is.numeric(lim) || is.na(lim) || lim < 1) lim <- 3L
  lim
}

.emc_stall_abort <- function(stage = NULL, iteration = NULL,
                             run_length = NULL, chains = NULL) {
  location <- paste0(
    if (!is.null(stage)) paste0(" stage=", stage) else "",
    if (!is.null(iteration)) paste0(" iteration=", iteration) else "",
    if (!is.null(chains)) paste0(" chain=", paste(chains, collapse = ",")) else ""
  )
  run_text <- if(is.null(run_length)) "" else {
    paste0(" (", run_length, " consecutive tries without moving)")
  }
  stop(errorCondition(
    sprintf(paste0("EMC2 aborted after a stalled chain%s%s.\n",
                   "  options(emc2.stall_try_limit = Inf) to disable the check."),
            location, run_text),
    class = c("emc_sampler_stalled", "emc_failure"),
    emc_class = "numerical", emc_source = "sampler",
    emc_stage = stage, emc_iteration = iteration,
    emc_run_length = run_length, emc_chains = chains))
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
