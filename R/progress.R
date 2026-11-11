# =============================================================================
# progress.R — the progress and cancel protocol (0.9.0, design D1 §3.5)
#
# One contract for every long call. The caller passes `progress = NULL` or a
# function of one argument; the package calls it with a plain event list
#
#   list(stage, unit, i, n, label, elapsed, payload)
#
# after each unit of work. The return value is ignored unless it is FALSE: then
# the package stops at that point with `hs_error_cancelled`, carrying
# `partial` (the work completed so far) and `at = list(unit, i, n)`.
#
# .hs_progress_tracker() is internal and never stored in a result. With
# `progress = NULL` a tick costs one is.null(). `partial` is a thunk, so nothing
# is built unless the run is cancelled.
# =============================================================================

.hs_progress_tracker <- function(progress, stage, unit, n) {
  if (!is.null(progress) && !is.function(progress))
    .hs_abort("'progress' must be NULL or a function of one argument (the event).",
              call. = FALSE)
  t0 <- Sys.time()
  n  <- as.integer(n)
  tick <- function(i, label = NA_character_, payload = NULL,
                   partial = function() NULL) {
    if (is.null(progress)) return(invisible(TRUE))
    ev <- list(stage = stage, unit = unit, i = as.integer(i), n = n,
               label = as.character(label)[1L],
               elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")),
               payload = payload)
    go <- progress(ev)
    if (isFALSE(go))
      .hs_abort(stage, ": cancelled by the progress callback after ", unit, " ",
                i, " of ", n, ".", call. = FALSE, class = "hs_error_cancelled",
                data = list(partial = partial(),
                            at = list(unit = unit, i = as.integer(i), n = n)))
    invisible(TRUE)
  }
  list(tick = tick, n = n, stage = stage, unit = unit)
}
