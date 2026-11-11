# =============================================================================
# state.R — process state the package touches, and how it gives it back
#
# A haploStack call must leave the R process as it found it (design framework
# F1/F2, 0.9.0): the RNG stream, the options, the working directory, the
# global environment, the future plan, the graphics devices. This file holds
# the shared helpers that make that true. Each is the only place its kind of
# state is touched:
#
#   .with_seed(seed, code)         evaluate `code` under set.seed(seed), then
#                                  put the caller's stream back
#   .local_seed(seed)              the same, scoped to the calling function
#                                  (for long bodies; registers the restore on
#                                  the caller's on.exit)
#   .with_rng_state(state, code)   evaluate under a captured stream and return
#                                  the new stream (resumable jobs)
#   .hs_with_plan(n_cores, code)   the future plan for a parallel section, put
#                                  back afterwards (D1-4)
#   .hs_keep_plan(code)            run `code` and put the caller's plan back
#                                  (a backend that replaces it, e.g. HapSelect)
#   .hs_keep_options(names, code)  run `code` and put the named options back
#
# No other file calls set.seed(), future::plan() or reads getOption() (hygiene
# tests, D1-3 / D1-4 / D1-5).
# =============================================================================

# Snapshot of the global RNG stream: NULL when the session has not drawn yet.
.hs_rng_get <- function() {
  genv <- globalenv()
  if (exists(".Random.seed", envir = genv, inherits = FALSE))
    get(".Random.seed", envir = genv, inherits = FALSE)
  else NULL
}

# Put a snapshot back. NULL restores the "not drawn yet" state, so the next
# draw seeds itself as it would have in an untouched session. `.Random.seed`
# carries the generator kind, so restoring it restores the kind too.
.hs_rng_set <- function(state) {
  genv <- globalenv()
  if (is.null(state)) {
    if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(list = ".Random.seed", envir = genv)
  } else {
    assign(".Random.seed", state, envir = genv)
  }
  invisible(NULL)
}

# Evaluate `code` with the RNG seeded to `seed`; the caller's stream is
# restored on exit, whether `code` returns or throws. `seed = NULL` evaluates
# on the caller's stream, untouched. A function's output for a given seed is
# exactly what a bare set.seed(seed) gave; only the state left behind differs.
.with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  old <- .hs_rng_get()
  on.exit(.hs_rng_set(old), add = TRUE)
  set.seed(seed)
  code
}

# `.with_seed()` for a whole function body: call it at the top of a function
# and the restore runs when that function exits. `seed = NULL` does nothing.
.local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  old <- .hs_rng_get()
  do.call(base::on.exit,
          list(substitute(.hs_rng_set(old), list(old = old)), TRUE, FALSE),
          envir = envir)
  set.seed(seed)
  invisible(NULL)
}

# Run `code` and put the caller's stream back afterwards, whatever `code` drew
# (an entry point whose jobs seed themselves, e.g. the simulation wrapper).
.hs_keep_rng <- function(code) {
  old <- .hs_rng_get()
  on.exit(.hs_rng_set(old), add = TRUE)
  code
}

# Evaluate `code` on a stream of its own and hand the stream back, so a job
# can be resumed later from exactly where it stopped. `state` is either a
# stream captured by an earlier call (the integer vector `.Random.seed`
# holds) or, for the first step, a seed. The caller's own stream is swapped
# out for the call and put back afterwards. Returns list(value, state).
.with_rng_state <- function(state, code) {
  old <- .hs_rng_get()
  on.exit(.hs_rng_set(old), add = TRUE)
  if (is.integer(state) && length(state) > 1L) .hs_rng_set(state)
  else set.seed(state)
  value <- code
  list(value = value, state = .hs_rng_get())
}


# =============================================================================
# The future plan (design D1 §3.7, build step D1-4)
#
# Three entry points run work in parallel when `n_cores > 1`
# (run_haplotype_catalogue(), run_stacking_simulation() twice). Up to 0.8.8 each
# set future::plan(multisession) itself and put the old plan back; a caller who
# had set a plan of their own (a cluster, a different worker count) had it
# overridden for the call. One helper now owns that decision.
# =============================================================================

# Can a parallel section run? TRUE when n_cores > 1 and both packages are
# installed; otherwise FALSE, with the 0.8.8 warning when parallel work was
# asked for and the packages are missing.
.hs_parallel_ok <- function(n_cores) {
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores <= 1L) return(FALSE)
  ok <- requireNamespace("future", quietly = TRUE) &&
        requireNamespace("future.apply", quietly = TRUE)
  if (!ok)
    .hs_warn("'future' and 'future.apply' required for n_cores > 1. ",
             "Falling back to sequential execution.", call. = FALSE)
  ok
}

# Evaluate `code` under the plan a parallel section needs, and put the caller's
# plan back on exit. `n_cores <= 1`: the plan is not touched. `n_cores > 1`
# with the caller's plan sequential (the default): multisession with `n_cores`
# workers for the call (0.8.8's behaviour, so scripts see no change). `n_cores
# > 1` with a plan the caller set themselves (anything not sequential): that
# plan is used as it is, and one message says so (agent decision, D1 §3.7).
.hs_with_plan <- function(n_cores, code) {
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores <= 1L) return(code)
  if (!requireNamespace("future", quietly = TRUE))
    .hs_abort(class = "hs_error_dependency",
              data = list(hint = "install.packages(\"future\")"),
              "n_cores > 1 needs the 'future' package.", call. = FALSE)
  current <- future::plan()
  if (inherits(current, "sequential")) {
    old <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(old), add = TRUE)
  } else {
    kind <- intersect(c("multisession", "multicore", "cluster", "sequential"),
                      class(current))[1L]
    if (is.na(kind)) kind <- class(current)[1L]
    .hs_inform("using the caller's future plan (", kind, ", ",
               future::nbrOfWorkers(), " workers) for n_cores = ", n_cores, ".")
  }
  code
}

# Run `code` and put the caller's future plan back whatever `code` did to it
# (HapSelect::pairwise_ld(parallelize = TRUE) replaces the plan and resets it).
.hs_keep_plan <- function(code) {
  if (requireNamespace("future", quietly = TRUE)) {
    old <- future::plan()
    on.exit(future::plan(old), add = TRUE)
  }
  code
}

# The current values of the named options, as a list options() can restore
# (an absent option comes back as NULL, which options() unsets again). The one
# place that reads options by name outside the settings resolver (hygiene H1).
.hs_options_get <- function(names) {
  old <- lapply(names, getOption)
  names(old) <- names
  old
}

# Run `code` and put the named options back as they were.
.hs_keep_options <- function(names, code) {
  old <- .hs_options_get(names)
  on.exit(options(old), add = TRUE)
  code
}


# =============================================================================
# Settings: argument first, option as default (framework F1, form row 2 a)
#
# Five process settings exist (`haploStack.*` options, see zzz.R). Each is
# resolved ONCE, at the exported entry point, through .hs_setting(), and then
# passed down as an ordinary argument; no helper reads getOption("haploStack.")
# itself (hygiene test, test-090-hygiene.R: this file and zzz.R are the only
# readers). Precedence: the call's own argument > the option > the built-in
# default. An argument is validated exactly as the setter validates it, so a
# bad value fails at the call the user can see.
#
#   name            option                      default  read by
#   missing_chars   haploStack.missing_chars    NULL     the build/repair portal
#                                                        only (.extra_missing_chars)
#   sample_order    haploStack.sample_order     NULL     .hs_sample_order()
#   fill_na         haploStack.fill_na          "-"      .hs_fill_na()
#   strict_missing  haploStack.strict_missing   FALSE    .notify_missing_tokens()
#   use_rcpp        haploStack.use_rcpp         TRUE     .hs_use_rcpp()
#
# use_rcpp stays a process option on purpose (D1 §3.2): it picks a code path,
# never a result, so it is not per-job state.
# =============================================================================

# Validators. Each takes the raw value, where it came from ("argument" or
# "option") and the setting name, and returns the value to use or raises
# hs_error_input. An option that is merely unusable (rather than wrong) is
# ignored where the 0.8.8 reader ignored it, so options() set by hand cannot
# change a result silently.
.check_flag <- function(x, source, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    .hs_abort("'", name, "' must be TRUE or FALSE (", source, ").", call. = FALSE)
  x
}

.check_fill_na <- function(x, source, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x))
    .hs_abort("fill_na must be a single non-NA character string; got ",
              paste(utils::capture.output(utils::str(x)), collapse = " "), ".",
              call. = FALSE)
  if (!nzchar(x))
    .hs_abort('fill_na must not be the empty string — an empty placeholder shortens ',
              'the signature and lets distinct allele combinations collide. Use "-" ',
              "(the default) or another single character.", call. = FALSE)
  x
}

# The registry option is a per-mode list (set by use_missing_chars()); a call
# may instead pass a character vector, meaning "for this call's mode". The
# characters themselves are checked against the mode in .extra_missing_chars().
.check_missing_reg <- function(x, source, name) {
  if (is.list(x)) {
    ok <- vapply(x, function(v) is.character(v) || is.null(v), logical(1L))
    if (length(x) && (is.null(names(x)) || any(!nzchar(names(x))) || !all(ok))) {
      if (identical(source, "option")) return(NULL)
      .hs_abort("'missing_chars' must be a character vector, or a list named by ",
                "mode (nucleotide / selectiontools / numeric) of character ",
                "vectors.", call. = FALSE)
    }
    return(x)
  }
  if (is.character(x)) {
    if (identical(source, "option")) return(NULL)   # the option is always a list
    return(x)
  }
  if (identical(source, "option")) return(NULL)
  .hs_abort("'missing_chars' must be a character vector, or a list named by ",
            "mode of character vectors; got ", class(x)[1L], ".", call. = FALSE)
}

# An unusable option is ignored with a warning (the 0.8.8 rule); an argument is
# validated by use_sample_order() itself.
.check_sample_order <- function(x, source, name) {
  if (inherits(x, "hs_sample_clustering")) return(x)
  if (is.character(x) || is.factor(x)) return(as.character(x))
  if (identical(source, "option")) {
    .hs_warn("haploStack.sample_order holds an unusable value (",
             class(x)[1L], "); ignoring it. Set it with use_sample_order().",
             call. = FALSE)
    return(NULL)
  }
  .hs_abort("'sample_order' must be an hs_sample_clustering (from ",
            "cluster_samples()), a character vector of sample IDs, or NULL.",
            call. = FALSE)
}

.HS_SETTINGS <- list(
  missing_chars  = list(opt = "haploStack.missing_chars",  default = NULL,
                        check = .check_missing_reg),
  sample_order   = list(opt = "haploStack.sample_order",   default = NULL,
                        check = .check_sample_order),
  fill_na        = list(opt = "haploStack.fill_na",        default = HS_FILL_NA_DEFAULT,
                        check = .check_fill_na),          # R/alphabet.R loads first
  strict_missing = list(opt = "haploStack.strict_missing", default = FALSE,
                        check = .check_flag),
  use_rcpp       = list(opt = "haploStack.use_rcpp",       default = TRUE,
                        check = .check_flag),
  # SelectionTools formats (R/st_formats.R): off until use_st_formats(TRUE)
  st_formats     = list(opt = "haploStack.st_formats",     default = FALSE,
                        check = .check_flag))

# The one resolver. `value` is the call's argument (NULL = not given).
.hs_setting <- function(name, value = NULL) {
  spec <- .HS_SETTINGS[[name]]
  if (is.null(spec))
    .hs_abort(class = "hs_error_internal", ".hs_setting: unknown setting '",
              name, "'.", call. = FALSE)
  if (!is.null(value)) return(spec$check(value, "argument", name))
  x <- getOption(spec$opt)
  if (is.null(x)) return(spec$default)
  x <- spec$check(x, "option", name)
  if (is.null(x)) spec$default else x
}

# The stored option as it is, unvalidated: what a setter returns as the
# previous value, so it can be put back with options().
.hs_option_raw <- function(name) {
  getOption(.HS_SETTINGS[[name]]$opt)
}


# ---- locale-independent ordering (D1 §3.11, decision-form row 33 a) ----------
# Identifiers (sample, marker, block, chromosome and allele names) are sorted
# in C-locale byte order wherever the package sorts them, so a Windows
# desktop, a Linux server and a Chinese-locale session order a table the same
# way. sort(method = "radix") / order(method = "radix") are that order (the one
# data.table already uses); numeric sorts are unaffected. Chromosome names go
# through hs_chrom_levels() (natural order) regardless.
.hs_sort <- function(x, decreasing = FALSE) {
  if (is.factor(x)) x <- as.character(x)
  sort(x, decreasing = decreasing, method = "radix")
}

.hs_order <- function(..., decreasing = FALSE) {
  order(..., decreasing = decreasing, method = "radix")
}


# ---- the server contract (D1 §3.1, decision-form row 1 a) -------------------

#' Running haploStack as a server backend
#'
#' The supported server model is \strong{one job per R worker process}
#' (\code{callr} per job, or a reused \code{mirai} / \code{callr} session).
#' Every exported compute function leaves the process as it found it, so a
#' worker can run job after job; a reused worker calls
#' \code{\link{hs_reset_session}} between jobs to be sure.
#'
#' @section Guarantees (0.9.0):
#' \tabular{lll}{
#'   \strong{#} \tab \strong{Guarantee} \tab \strong{Mechanism} \cr
#'   G1 \tab No exported compute function leaves process state changed: the
#'     random stream, the options, the working directory, the global
#'     environment, the future plan, the graphics devices, circlize's
#'     parameters, files under \code{tempdir()}, an open genomicSimulation
#'     session \tab local seeds (\code{.with_seed()}), the settings resolver,
#'     \code{.hs_with_plan()}, \code{\link{hs_render_figure}}, private
#'     scratch removed on exit, the SelectionTools io-dir shim; the leak test
#'     \code{test-090-no-leak.R} \cr
#'   G2 \tab Every setting that changes a result can be passed as an
#'     argument; a worker never needs the setters \tab \code{missing_chars =},
#'     \code{center =}, \code{fill_na =}, \code{strict_missing =},
#'     \code{st_formats =} on the entry points \cr
#'   G3 \tab Errors are classed so an API can map them; a cancelled job
#'     returns its partial result \tab \code{hs_error_input} /
#'     \code{_format} / \code{_schema} / \code{_dependency} / \code{_state}
#'     / \code{_cancelled} (\code{e$partial}); the \code{progress} callback
#'     (\code{FALSE} = cancel) \cr
#'   G4 \tab genomicSimulation: one open session per process, a second load
#'     refused; parallel replications run in their own processes \tab
#'     \code{\link{gs_load_geno_meta}}, \code{\link{stk_init}} \cr
#'   G5 \tab User uploads go through the text readers; an \code{.rds} is
#'     read only from the site's own store \tab \code{\link{load_geno_meta}}
#'     (trusted input only) \cr
#'   G6 \tab Output ordering does not depend on the worker's locale \tab
#'     C-locale ordering of identifiers, \code{\link{hs_chrom_levels}} for
#'     chromosomes \cr
#' }
#'
#' @section Backends on a server:
#' HapSelect and genomicSimulation (GitHub, compiled) are supported in a
#' worker; SelectionTools stays off the server (decision-form row 31 a): its
#' thread count, resident LD data set and global variables cannot be released
#' by the package. rrBLUP, BGLR, circlize, cowplot and ragg (headless PNG) are
#' the CRAN extras a worker image needs; shiny and zip are not.
#'
#' @name haploStack_server
#' @seealso \code{\link{hs_reset_session}}, \code{\link{hs_quiet}},
#'   \code{\link{stk_step}}
NULL

#' Reset the package's process state in a reused worker
#'
#' Puts back what a haploStack call could have left behind in a reused R
#' process (see \link{haploStack_server}): clears the per-job settings
#' (\code{haploStack.sample_order}, \code{haploStack.missing_chars},
#' \code{haploStack.strict_missing}), resets the switches
#' (\code{haploStack.use_rcpp}, \code{haploStack.fill_na},
#' \code{haploStack.st_formats}) to their defaults, closes an open
#' genomicSimulation session, clears circlize's parameters when circlize is
#' loaded, and removes the package's own scratch directories under
#' \code{tempdir()} (prefix \code{hs_}). It never touches graphics devices,
#' the random stream, the future plan or files it did not create. With one
#' fresh process per job it is not needed.
#'
#' @param quiet When \code{FALSE}, a message lists what was reset.
#' @return Invisibly, a list of what was reset: \code{options} (the option
#'   names), \code{gs_session} (\code{TRUE} when one was closed),
#'   \code{circlize} (\code{TRUE} when cleared), \code{scratch} (the
#'   directory names removed).
#' @examples
#' \dontrun{
#' hs_reset_session()
#' }
#' @export
hs_reset_session <- function(quiet = TRUE) {
  done <- list(options = character(0L), gs_session = FALSE, circlize = FALSE,
               scratch = character(0L))
  # settings: the per-job ones are cleared, the switches go back to defaults
  clear <- c("missing_chars", "sample_order", "strict_missing")
  reset <- c("fill_na", "use_rcpp", "st_formats")
  opts  <- list()
  for (nm in clear) {
    o <- .HS_SETTINGS[[nm]]$opt
    if (!is.null(getOption(o))) { opts[o] <- list(NULL); done$options <- c(done$options, o) }
  }
  for (nm in reset) {
    o <- .HS_SETTINGS[[nm]]$opt
    d <- .HS_SETTINGS[[nm]]$default
    if (!identical(getOption(o), d)) { opts[[o]] <- d; done$options <- c(done$options, o) }
  }
  if (length(opts)) options(opts)
  # genomicSimulation: one session per process; close it
  if (isTRUE(.hs_gs$open)) { .gs_session_clear(); done$gs_session <- TRUE }
  # circlize keeps its parameters in its namespace until circos.clear()
  if ("circlize" %in% loadedNamespaces()) {
    try(circlize::circos.clear(), silent = TRUE)
    done$circlize <- TRUE
  }
  # the package's own scratch (tempfile("hs_<what>_") everywhere)
  dirs <- list.files(tempdir(), pattern = "^hs_", full.names = TRUE)
  if (length(dirs)) {
    unlink(dirs, recursive = TRUE, force = TRUE)
    done$scratch <- basename(dirs)
  }
  if (!isTRUE(quiet))
    .hs_inform("hs_reset_session: ",
               if (length(done$options)) paste0("options ", paste(done$options, collapse = ", "), "; ")
               else "no option to reset; ",
               if (done$gs_session) "genomicSimulation session closed; " else "",
               if (done$circlize) "circlize cleared; " else "",
               length(done$scratch), " scratch director",
               if (length(done$scratch) == 1L) "y" else "ies", " removed.")
  invisible(done)
}
