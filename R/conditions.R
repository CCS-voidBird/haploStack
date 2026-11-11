# =============================================================================
# conditions.R — classed errors, warnings and messages (0.9.0)
#
# Every error, warning and message the package raises goes through one of the
# three helpers below. They build the same text and the same call as
# base::stop() / warning() / message() did up to 0.8.8, and add a class:
#
#   error   c(<class>, "hs_error",   "error",   "condition")
#   warning c(<class>, "hs_warning", "warning", "condition")
#   message c("hs_message", "message", "condition")
#
# so a caller — a website worker in particular — can tell a haploStack
# condition from an R-level failure (`tryCatch(..., hs_error = ...)`), map
# input errors to a client response, and silence the package's progress chatter
# without silencing anything else (hs_quiet()).
#
# Class vocabulary (design framework F3 / A3). Every error is an input error
# unless it names another kind; input subkinds are written as
# c("hs_error_block_map", "hs_error_input") so both classes match.
#
#   hs_error_input        the caller's data or arguments are wrong   (default)
#   hs_error_block_map    ⊂ input: malformed block map (chromosome mismatch,
#                         overlap, no marker matched)
#   hs_error_format       ⊂ input: unreadable file, undetectable encoding,
#                         bad header
#   hs_error_schema       ⊂ input: container from a newer schema, or not
#                         upgradable
#   hs_error_defunct      ⊂ input (also base "defunctError"): a removed name
#   hs_error_dependency   a Suggests package needed for this call is absent;
#                         carries `hint` = the install command
#   hs_error_state        process state forbids the call (second
#                         genomicSimulation session, SelectionTools globals)
#   hs_error_cancelled    a progress callback asked to stop; carries `partial`
#                         (completed work) and `at` (unit, i, n)
#   hs_error_internal     an assertion failed; a package bug, not the caller's
#   hs_warning_input      a data problem that does not stop the call
#   hs_warning_fallback   a silent fallback that now leaves a trace (also a
#                         row in a result's notes); carries `stage`, `code`
#   hs_warning_deprecated deprecated name or argument (also base
#                         "deprecatedWarning")
#   hs_message            progress chatter; muffled by hs_quiet()
#
# Structured fields: `.hs_abort()` / `.hs_warn()` take `data = list(...)`,
# merged into the condition object (reserved names: hint, partial, at, bad).
# The message text and the call are never changed by `data`.
# =============================================================================

# Build a condition object with extra fields. `ctor` is errorCondition or
# warningCondition. `data` must be a named list; the names "message", "call"
# and "class" belong to the condition itself.
.hs_condition <- function(msg, call, class, data, ctor) {
  if (length(data)) {
    nm <- names(data)
    if (is.null(nm) || any(!nzchar(nm)))
      stop(errorCondition(".hs_condition: every 'data' field must be named.",
                          class = c("hs_error_internal", "hs_error"), call = NULL))
    if (any(nm %in% c("message", "call", "class")))
      stop(errorCondition(paste0(".hs_condition: 'data' may not use the reserved ",
                                 "field names message, call or class."),
                          class = c("hs_error_internal", "hs_error"), call = NULL))
  }
  # quote = TRUE: `call` is a language object, and do.call() would otherwise
  # evaluate it (re-running the failing call) while building the condition.
  do.call(ctor, c(list(msg, class = class, call = call), data), quote = TRUE)
}

.hs_abort <- function(..., call. = TRUE, domain = NULL, class = "hs_error_input",
                      data = list()) {
  if (...length() == 1L && inherits(..1, "condition")) stop(..1)
  msg  <- .makeMessage(..., domain = domain)
  call <- if (isTRUE(call.)) sys.call(-1L) else NULL
  stop(.hs_condition(msg, call, c(class, "hs_error"), data, errorCondition))
}

.hs_warn <- function(..., call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
                     domain = NULL, class = character(), data = list()) {
  if (...length() == 1L && inherits(..1, "condition")) return(warning(..1))
  msg  <- .makeMessage(..., domain = domain)
  call <- if (isTRUE(call.)) sys.call(-1L) else NULL
  if (isTRUE(immediate.) || isTRUE(noBreaks.)) {
    # Classed warnings cannot carry these flags; keep the base behaviour.
    return(warning(msg, call. = FALSE, immediate. = immediate.,
                   noBreaks. = noBreaks.))
  }
  warning(.hs_condition(msg, call, c(class, "hs_warning"), data, warningCondition))
}

.hs_inform <- function(..., domain = NULL, appendLF = TRUE) {
  if (...length() == 1L && inherits(..1, "condition")) return(message(..1))
  msg  <- .makeMessage(..., domain = domain, appendLF = appendLF)
  call <- sys.call(-1L)
  message(structure(class = c("hs_message", "message", "condition"),
                    list(message = msg, call = call)))
}


#' Evaluate an expression without haploStack's progress messages
#'
#' Muffles the informational messages haploStack raises (class
#' \code{"hs_message"}) while letting every other message, and every warning
#' and error, through. Useful in scripts and in a server worker, where the
#' package's progress lines would otherwise flood the log. Nothing global is
#' changed: the effect ends with the call.
#'
#' Since 0.9.0 the package's conditions carry classes. Every error is
#' \code{"hs_error"}; by default it is also \code{"hs_error_input"} (the
#' caller's data or arguments are wrong), and other kinds replace that:
#' \code{"hs_error_dependency"} (an optional package is missing; the condition's
#' \code{hint} field holds the install command), \code{"hs_error_state"} (the
#' process state forbids the call), \code{"hs_error_cancelled"} (a progress
#' callback stopped the call; \code{partial} holds the completed work) and
#' \code{"hs_error_internal"} (a package assertion failed). Input errors may
#' carry a finer class first, such as \code{"hs_error_block_map"},
#' \code{"hs_error_format"} or \code{"hs_error_schema"}. Warnings are
#' \code{"hs_warning"} (with \code{"hs_warning_input"},
#' \code{"hs_warning_fallback"} or \code{"hs_warning_deprecated"} in front),
#' messages \code{"hs_message"}. Use them with \code{tryCatch()} or
#' \code{withCallingHandlers()} to treat haploStack's conditions apart from
#' R's own.
#'
#' @param expr Expression to evaluate.
#' @return The value of \code{expr}.
#' @examples
#' \dontrun{
#' meta <- hs_quiet(read_plink_ped("panel.ped", "panel.map"))
#' res  <- tryCatch(run_haplotype_catalogue(meta, output_dir = NULL),
#'                  hs_error = function(e) conditionMessage(e))
#' }
#' @export
hs_quiet <- function(expr) {
  withCallingHandlers(expr, hs_message = function(m) invokeRestart("muffleMessage"))
}
