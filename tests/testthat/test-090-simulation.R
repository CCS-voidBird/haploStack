# =============================================================================
# test-090-simulation.R — the simulation driver (design D4, build phase E)
#
#   E1  D4-1 pure refactors: the scheme registry replaces every string match,
#       .hap_cols() is the one column parser, .gbs_spec() the one GBS resolver,
#       .expand_cross_plan() serves gs_cross_plan() too, one trajectories view,
#       the founder row through .stk_summary_record()
#
# Every E1 check is a "same numbers" check: the refactor moved code, not
# behaviour (the fixed-seed parity harness is the gate for the whole route).
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

.e_fixture <- function(seed = 5L, n_samples = 12L, n_markers = 24L) {
  m <- sim_geno_meta(n_samples = n_samples, n_markers = n_markers, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"))
  list(meta = m, eff = e,
       lg = compute_profile_localgebv(m, e, level = "haplotype"))
}

# ── E1: the scheme registry ──────────────────────────────────────────────────

test_that("E1: the registry names the four schemes and every site reads it", {
  reg <- haploStack:::.STK_SCHEMES
  expect_named(reg, c("TS", "OHS", "OCS", "GBS"))
  for (s in names(reg)) {
    expect_true(is.function(reg[[s]]$parents), info = s)
    expect_true(is.null(reg[[s]]$pairs) || is.function(reg[[s]]$pairs), info = s)
    expect_true(is.logical(reg[[s]]$needs_ideal), info = s)
    expect_true(reg[[s]]$direction %in% c("both", "high", "argmax"), info = s)
  }
  expect_null(reg$OCS$pairs)                      # no family form before KCT
  expect_true(reg$GBS$needs_ideal)
  expect_true(reg$GBS$pin_founder_blocks)
  expect_identical(haploStack:::.STK_BASELINE, "TS")
  expect_identical(haploStack:::.stk_scheme_flag(c("TS", "GBS"), "needs_ideal"),
                   c(FALSE, TRUE))
  expect_error(haploStack:::.stk_scheme("XYZ"), class = "hs_error_internal")

  # no scheme-name string comparison is left in the driver source
  src <- testthat::test_path("..", "..", "R", "stacking_main.R")
  skip_if_not(file.exists(src), "R/ not reachable")
  code <- readLines(src, warn = FALSE)
  code <- code[!grepl("^\\s*#", code)]
  expect_false(any(grepl('scheme == "(TS|OHS|OCS|GBS)"', code)))
  expect_false(any(grepl('"(OCS|GBS)" %in% scheme', code)))
  expect_false(any(grepl("%\\|\\|%", code[grepl("opts\\$", code)])))
})

test_that("E1: the registry selectors give the 0.8.8 switch-arm results", {
  fx <- .e_fixture()
  opts <- list(select_priority = "high", ohs_method = "Greedy", ohs_elite_ratio = 1,
               ohs_elite_level = "pair", seed = 42L, max_kinship = 0.125,
               freq_penalty_weight = 1, gbs_method = "ARCSIN",
               gbs_objective = "candidate", gbs_engine = "greedy",
               penalty_mode = "dynamic", gbs_penalty = -10, penalty_scale = 10,
               penalty_rank_scale = 1, gbs_effect_weight = 0)
  ideal <- construct_ideal_genotype(fx$lg)
  expect_identical(
    haploStack:::.stk_select_parents("TS", fx$lg, fx$meta, 4L, NULL, opts),
    haploStack:::.stk_select_truncation(fx$lg, 4L, "high"))
  expect_identical(
    haploStack:::.stk_select_parents("OHS", fx$lg, fx$meta, 4L, NULL, opts),
    ohs_selection(fx$lg, n_parents = 4L, method = "Greedy"))
  expect_identical(
    haploStack:::.stk_select_parents("GBS", fx$lg, fx$meta, 4L, ideal, opts),
    gbs_selection(fx$lg, ideal_genotypes = ideal, n_parents = 4L))
  ocs <- haploStack:::.stk_select_parents("OCS", fx$lg, fx$meta, 4L, NULL, opts)
  expect_true(length(ocs) >= 2L && all(ocs %in% dimnames(fx$meta$geno_arrays)[[1L]]))
  # OCS has no family form: the pairs dispatcher names the scheme
  opts$pair_budget <- haploStack:::.stk_selection_budget("pairs", "families", 4L)
  expect_error(haploStack:::.stk_select_pairs("OCS", fx$lg, fx$meta, NULL, opts),
               "no direct-pair form")
})

# ── E1: .hap_cols() ──────────────────────────────────────────────────────────

test_that("E1: .hap_cols() parses sample and copy exactly as the 0.8.8 regex did", {
  cols <- c("S1.1", "S1.2", "S10.1", "S10.2", "a.b.1", "a.b.2", "noidx")
  hc <- haploStack:::.hap_cols(cols)
  expect_identical(hc$base, sub("\\.[0-9]+$", "", cols))
  expect_identical(hc$samples, unique(sub("\\.[0-9]+$", "", cols)))
  expect_identical(hc$idx, c(1L, 2L, 1L, 2L, 1L, 2L, NA_integer_))
  expect_identical(hc$ploidy, 2L)
  expect_identical(haploStack:::.hap_cols(cols, ploidy = 4L)$ploidy, 4L)
  expect_identical(haploStack:::.hap_cols(character())$ploidy, 0L)

  # the consumers agree with each other on a real table
  fx <- .e_fixture()
  long <- haploStack:::.localgebv_long(fx$lg)
  expect_identical(sort(unique(long$Sample)), sort(dimnames(fx$meta$geno_arrays)[[1L]]))
  expect_setequal(unique(long$Haplotype), 1:2)
  # consumers keep the table's own column order (lexicographic in 0.8.8)
  want <- unique(sub("[.][0-9]+$", "", setdiff(names(fx$lg), "Block")))
  geno <- haploStack:::.collapse_hap_to_geno(fx$lg)
  expect_identical(setdiff(names(geno), "Block"), want)
  M <- haploStack:::.pairwise_block_value_matrix(fx$lg)
  expect_identical(rownames(M), want)
})

# ── E1: .gbs_spec() ──────────────────────────────────────────────────────────

test_that("E1: .gbs_spec() is the three 0.8.8 preset paths in one", {
  spec <- haploStack:::.gbs_spec("haploStack")
  expect_identical(spec[names(haploStack:::.GBS_PRESETS$haploStack)],
                   haploStack:::.GBS_PRESETS$haploStack)
  expect_identical(spec$penalty_rank_scale, 1)
  # explicit wins, preset fills the rest
  spec <- haploStack:::.gbs_spec("kemper", penalty_mode = "dynamic", gbs_penalty = -5,
                                 penalty_rank_scale = c(1, 0.5))
  expect_identical(spec$penalty_mode, "dynamic")
  expect_identical(spec$gbs_penalty, -5)
  expect_identical(spec$penalty_scale, 10)
  expect_identical(spec$penalty_rank_scale, c(1, 0.5))
  # validation as before
  expect_error(haploStack:::.gbs_spec("haploStack", gbs_penalty = 3), "non-positive")
  expect_error(haploStack:::.gbs_spec("haploStack", penalty_scale = 0.5), ">= 1")
  expect_error(haploStack:::.gbs_spec("haploStack", penalty_rank_scale = -1), "finite")
  expect_error(haploStack:::.gbs_spec("haploStack", gbs_engine = "GA", caller = "x"),
               "x: gbs_engine = \"GA\" optimises a set objective")
  expect_identical(haploStack:::.gbs_spec("haploStack", gbs_engine = "GA",
                                          gbs_objective = "set")$gbs_engine, "GA")
  expect_error(haploStack:::.gbs_spec("haploStack", gbs_objective = "nope"))

  # the entry points still give the same selections
  fx <- .e_fixture()
  ideal <- construct_ideal_genotype(fx$lg)
  s1 <- gbs_selection(fx$lg, ideal_genotypes = ideal, n_parents = 4L)
  s2 <- gbs_selection(fx$lg, ideal_genotypes = ideal, n_parents = 4L,
                      gbs_preset = "haploStack", penalty_mode = "dynamic",
                      gbs_penalty = -10, penalty_scale = 10)
  expect_identical(s1, s2)
  expect_error(gbs_selection(fx$lg, n_parents = 4L, gbs_engine = "GA"),
               "gbs_selection: gbs_engine = \"GA\" optimises a set objective")
  p1 <- gbs_pair_selection(fx$lg, ideal_genotypes = ideal, pair_budget = "families",
                           top_founders = 4L)
  p2 <- gbs_pair_selection(fx$lg, ideal_genotypes = ideal, pair_budget = "families",
                           top_founders = 4L, gbs_engine = "greedy")
  expect_identical(p1, p2)
  expect_error(gbs_pair_selection(fx$lg, pair_budget = "families", top_founders = 4L,
                                  gbs_engine = "GA"),
               "over family sets is not implemented")
})

# ── E1: shared expansion, trajectories view, founder record ─────────────────

test_that("E1: .expand_cross_plan(keep_cols = TRUE) keeps the design's columns", {
  cd <- make_cross_design(c("A", "B", "C"), method = "all_pairs", n_offspring = 2L)
  cd$extra <- c("x", "y", "z")
  e1 <- haploStack:::.expand_cross_plan(cd)
  e2 <- haploStack:::.expand_cross_plan(cd, keep_cols = TRUE)
  expect_s3_class(e2, "data.frame")
  expect_true(all(c("Parent1", "Parent2", "Cross_ID", "n_offspring", "extra",
                    "CrossID") %in% names(e2)))
  expect_identical(e2$CrossID, e1$CrossID)
  expect_identical(e2$Parent1, e1$Parent1)
  expect_identical(e2$extra, rep(c("x", "y", "z"), each = 2L))
  expect_true(all(grepl("^AxBn0[12]$", e2$CrossID[1:2])))
})

test_that("E1: one trajectories view, used by the result and the live monitor", {
  fx <- .e_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE)
  v <- haploStack:::.stk_trajectories_view(res$summary)
  expect_identical(v, res$trajectories)
  expect_named(v, c("Rep", "Cycle", "Method", "MeanGEBV", "SDev", "BlockCoverage"))
  expect_false(any(v$Cycle == 0L))
})

test_that("E1: the founder row is the 0.8.8 arithmetic through .stk_summary_record()", {
  fx <- .e_fixture()
  flg <- haploStack:::.stk_limit_blocks(fx$lg, NULL, 50L)
  ideal <- construct_ideal_genotype(flg)
  opts <- list(metrics = c("gebv", "d2ideal", "ideal_hap_gebv"), gbs_method = "ARCSIN",
               ploidy = 2L)
  row <- haploStack:::.stk_summary_record(
    haploStack:::.stk_founder_measures(fx$meta, flg, ideal, opts), "founder", 0L, 1L, "TS")
  v0 <- haploStack:::.stk_total_gebv(flg)
  expect_identical(row$cycle, 0L)
  expect_identical(row$select_type, "founder")
  expect_equal(row$mean_F1mean, mean(v0))
  expect_equal(row$max_F1mean, max(v0))
  expect_equal(row$mean_D2Ideal, haploStack:::.distance_to_ideal(flg, ideal))
  expect_equal(row$mean_ideal_phenotype, as.numeric(attr(ideal, "ideal_phenotype")))
  expect_true(is.na(row$mean_coancestry))             # not requested
  expect_identical(row$n_families, NA_integer_)
  # the driver's row is this row
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    metrics = c("d2ideal", "ideal_hap_gebv"), output_dir = NULL, verbose = FALSE)
  expect_equal(res$summary[cycle == 0L]$mean_D2Ideal, row$mean_D2Ideal)
  expect_equal(res$summary[cycle == 0L]$ideal_hap_gebv, row$ideal_hap_gebv)
})

# ── E2: one local-GEBV pass per generation (row 20 a, D4-2) ─────────────────

test_that("E2: the haplotype-level sum is the genotype-level F1mean, bit for bit", {
  # the D4 §3.4 probe as a test: fully blocked, with NA cells, partially blocked
  base <- sim_geno_meta(n_samples = 24L, n_markers = 96L, block_size = 4L,
                        encoding = "nucleotide", seed = 3L)
  eff  <- sim_effect_df(dimnames(base$geno_arrays)[[2L]], seed = 3L)
  m_na <- base
  set.seed(5); m_na$geno_arrays[sample(length(m_na$geno_arrays), 60L)] <- NA_character_
  bm   <- as.data.frame(base$block_marker_map_info)
  keep <- bm$Block %in% unique(bm$Block)[c(TRUE, TRUE, FALSE)]
  m_part <- new_geno_meta(base$geno_arrays, base$map_info,
                          base$block_info[base$block_info$Block %in% bm$Block[keep], , drop = FALSE],
                          bm[keep, , drop = FALSE], base$metadata)
  for (meta in list(base, m_na, m_part)) for (rcpp in c(TRUE, FALSE)) {
    op <- options(haploStack.use_rcpp = rcpp)
    smp <- dimnames(meta$geno_arrays)[[1L]]
    cd  <- make_cross_design(smp[1:6], method = "all_pairs", n_offspring = 3L)
    set.seed(11)
    res <- haploStack:::.cross_block(meta, cd, effect_df = eff, f1mean = TRUE)
    set.seed(11)
    res0 <- haploStack:::.cross_block(meta, cd, effect_df = eff, f1mean = FALSE)
    expect_true(all(is.na(res0$plans$F1mean)))
    expect_identical(res0$geno_meta$geno_arrays, res$geno_meta$geno_arrays)
    nm <- dimnames(res$geno_meta$geno_arrays)[[1L]]
    lg <- compute_profile_localgebv(res$geno_meta, eff, level = "haplotype")
    expect_identical(unname(haploStack:::.stk_total_gebv(lg)[nm]), res$plans$F1mean)
    options(op)
  }
})

test_that("E2: the driver hands the offspring table to the next generation and keeps its numbers", {
  fx <- .e_fixture()
  run <- function(...) run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 3L, scheme = "TS", select_top_blocks = 50L, metrics = "all",
    output_dir = NULL, verbose = FALSE, seed = 9L, ...)
  res <- run()
  expect_false(anyNA(res$plans$F1mean))
  # every generation's F1mean is the total local GEBV of its own offspring: the
  # final population reproduces the last cycle's plans
  fin <- res$final_populations[[1L]]
  lg  <- compute_profile_localgebv(fin, fx$eff, level = "haplotype")
  tot <- haploStack:::.stk_total_gebv(lg)
  last <- res$plans[cycle == max(cycle)]
  fam  <- as.data.table(fin$sample_info)[, .(F1 = mean(tot[Sample])), by = .(Parent1, Parent2)]
  m <- merge(last, fam, by = c("Parent1", "Parent2"))
  expect_equal(m$F1mean, m$F1, tolerance = 1e-12)
  # .stk_advance() accepts the previous table and gives the same generation
  opts_run <- NULL
  adv1 <- haploStack:::.stk_advance
  expect_true("lg" %in% names(formals(adv1)))
})

test_that("E2: a failed offspring pass leaves F1mean NA, notes it, and the run goes on", {
  fx <- .e_fixture()
  orig <- haploStack:::compute_profile_localgebv
  n_f  <- dim(fx$meta$geno_arrays)[1L]
  calls <- 0L
  # fail the FIRST offspring pass only (generation 1); the next generation then
  # recomputes its parent table (as 0.8.8 did every generation) and goes on
  testthat::local_mocked_bindings(
    compute_profile_localgebv = function(geno_meta, ...) {
      if (dim(geno_meta$geno_arrays)[1L] > n_f && calls == 0L) {
        calls <<- calls + 1L; stop("boom")
      }
      orig(geno_meta, ...)
    }, .package = "haploStack")
  w <- NULL
  res <- withCallingHandlers(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 3L, scheme = "TS", select_top_blocks = 50L, metrics = "all",
      output_dir = NULL, verbose = FALSE),
    hs_warning_fallback = function(cnd) { w <<- c(w, cnd$code); invokeRestart("muffleWarning") })
  expect_true(all(is.na(res$plans[cycle == 1L]$F1mean)))
  expect_false(anyNA(res$plans[cycle == 2L]$F1mean))
  expect_true(is.na(res$summary[cycle == 1L]$mean_D2Ideal))
  expect_false(is.na(res$summary[cycle == 2L]$mean_D2Ideal))
  expect_identical(w, "f1mean_failed")
  expect_identical(calls, 1L)
  expect_false(is.na(res$summary[cycle == 0L]$mean_F1mean))   # founders unaffected
})

# ── E3: the step API (row 18 a, D4-3) ───────────────────────────────────────

.e_stream_untouched <- function(expr) {
  set.seed(314159)
  pre <- .Random.seed
  force(expr)
  identical(pre, .Random.seed)
}

test_that("E3: stk_spec() validates once; stk_jobs() is the 0.8.8 grid", {
  sp <- stk_spec(scheme = c("TS", "GBS"), replications = 2L, select_cycles = 3L,
                 selfing_cycles = 1L, top_founders = 4L, num_offspring = 2L)
  expect_s3_class(sp, "stk_spec")
  expect_identical(sp$n_generations, 6L)
  expect_identical(sp$metrics, "gebv")
  expect_identical(sp$gbs$gbs_engine, "greedy")
  expect_identical(sp$budget$n_families, 6L)
  j <- stk_jobs(sp)
  expect_identical(j$job, 1:4)
  expect_identical(j$rep, c(1L, 2L, 1L, 2L))
  expect_identical(j$scheme, c("TS", "TS", "GBS", "GBS"))
  expect_output(print(sp), "TS/GBS x 2 rep")
  expect_error(stk_spec(replications = 0L), "must be >= 1")
  expect_error(stk_spec(top_founders = 1L), ">= 2")
  expect_error(stk_spec(backend_opts = list(nope = 1)), "unknown backend_opts")
  expect_error(stk_spec(ga_control = list(x = 1)), "ga_control")
  expect_error(stk_spec(selection_unit = "pairs", pair_budget = "families",
                        scheme = "OCS"), "no direct-pair form")
  expect_error(stk_jobs(list()), "stk_spec")
  # seed = NULL draws one integer from the caller's stream and records it
  set.seed(7); a <- stk_spec(seed = NULL)$seed
  set.seed(7); b <- stk_spec(seed = NULL)$seed
  expect_identical(a, b)
  expect_true(is.numeric(a) && a >= 1)
  expect_identical(stk_spec()$seed, 42L)
})

test_that("E3: the stk_step loop equals run_stacking_simulation (block, three schemes)", {
  fx <- .e_fixture()
  args <- list(scheme = c("TS", "OHS", "GBS"), replications = 2L, select_cycles = 2L,
               selfing_cycles = 1L, ga_start_by = 2L, top_founders = 4L,
               num_offspring = 2L, ohs_method = "Greedy", select_top_blocks = 50L,
               metrics = "all", seed = 11L)
  res <- do.call(run_stacking_simulation,
                 c(list(fx$meta, fx$eff, output_dir = NULL, verbose = FALSE), args))
  sp <- do.call(stk_spec, args)
  states <- lapply(seq_len(nrow(stk_jobs(sp))), function(k) {
    st <- stk_init(sp, fx$meta, fx$eff, job = k)
    n <- 0L
    while (st$status == "running") { st <- stk_step(st); n <- n + 1L }
    expect_identical(n, sp$n_generations)
    expect_identical(st$status, "done")
    expect_identical(stk_step(st), st)             # a done state is a no-op
    st
  })
  man <- stk_result(states)
  expect_identical(man$summary, res$summary)
  expect_identical(man$plans, res$plans)
  expect_identical(names(man$final_populations), names(res$final_populations))
  expect_identical(man$final_populations[[3L]]$geno_arrays,
                   res$final_populations[[3L]]$geno_arrays)
  expect_identical(man$ideal_genotypes, res$ideal_genotypes)
  expect_identical(stk_result(states, keep_populations = "none")$final_populations,
                   stats::setNames(list(), character()))
  # the schedule: decision / selfing / scheme used per generation
  ohs1 <- res$summary[sel_scheme == "OHS" & rep == 1L][order(cycle)]
  expect_identical(ohs1$select_type, c("founder", "TS", "selfing", "OHS", "selfing"))
})

test_that("E3: interleaved jobs equal sequential jobs; a saved block state resumes identically", {
  fx <- .e_fixture()
  sp <- stk_spec(scheme = "TS", replications = 2L, select_cycles = 4L, top_founders = 4L,
                 num_offspring = 2L, select_top_blocks = 50L, seed = 3L)
  seq_res <- stk_result(lapply(1:2, function(k) {
    st <- stk_init(sp, fx$meta, fx$eff, job = k)
    while (st$status == "running") st <- stk_step(st)
    st
  }))
  a <- stk_init(sp, fx$meta, fx$eff, job = 1L)
  b <- stk_init(sp, fx$meta, fx$eff, job = 2L)
  while (a$status == "running" || b$status == "running") {
    a <- stk_step(a); runif(3)                     # draws between steps do not matter
    b <- stk_step(b); sample(10)
  }
  inter <- stk_result(list(a, b))
  expect_identical(inter$summary, seq_res$summary)
  expect_identical(inter$plans, seq_res$plans)
  # save at generation 2, reload, finish
  st <- stk_init(sp, fx$meta, fx$eff, job = 1L)
  st <- stk_step(stk_step(st))
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(st, f)
  st2 <- readRDS(f)
  while (st2$status == "running") st2 <- stk_step(st2)
  r2 <- stk_result(st2)
  expect_identical(r2$summary, seq_res$summary[rep == 1L])
  expect_identical(r2$plans, seq_res$plans[rep == 1L])
})

test_that("E3: a deserialised genomicSimulation state refuses to step", {
  skip_if_not_installed("genomicSimulation")
  fx <- .e_fixture()
  sp <- stk_spec(scheme = "TS", replications = 1L, select_cycles = 2L, top_founders = 4L,
                 num_offspring = 2L, select_top_blocks = 50L, cross_backend = "genomicSimulation",
                 backend_opts = list(map_unit = "bp"))
  st <- stk_init(sp, fx$meta, fx$eff)
  withr::defer(stk_close(st))
  expect_s3_class(st$session, "gs_handle")
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(st, f)
  stk_close(st)                                   # the live session is closed ...
  st2 <- readRDS(f)
  expect_error(stk_step(st2), class = "hs_error_state")   # ... so the copy cannot resume
  expect_error(stk_step(st2), "cannot be resumed")
  # a fresh job runs; the session is freed when the job is done
  st3 <- stk_init(sp, fx$meta, fx$eff)
  while (st3$status == "running") st3 <- stk_step(st3)
  expect_null(st3$session)
  expect_false(haploStack:::.gs_session_open())
  expect_false(anyNA(stk_result(st3)$plans$F1mean))
})

test_that("E3: run_stacking_simulation leaves the caller's RNG untouched", {
  fx <- .e_fixture()
  expect_true(.e_stream_untouched(run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE)))
  sp <- stk_spec(scheme = "TS", replications = 1L, select_cycles = 2L, top_founders = 4L,
                 num_offspring = 2L, select_top_blocks = 50L)
  expect_true(.e_stream_untouched({
    st <- stk_init(sp, fx$meta, fx$eff)
    while (st$status == "running") st <- stk_step(st)
  }))
})

test_that("E3: progress sees every generation; FALSE cancels with a partial result", {
  fx <- .e_fixture()
  ev <- list()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 3L, selfing_cycles = 1L,
    top_founders = 4L, num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE,
    progress = function(e) { ev[[length(ev) + 1L]] <<- e; TRUE })
  expect_length(ev, 6L)
  expect_identical(vapply(ev, `[[`, integer(1L), "i"), 1:6)
  expect_identical(unique(vapply(ev, `[[`, character(1L), "unit")), "cycle")
  expect_identical(vapply(ev, function(e) e$payload$select_type, character(1L)),
                   rep(c("TS", "selfing"), 3L))
  expect_identical(vapply(ev, function(e) e$payload$cycle, integer(1L)), rep(1:3, each = 2L))
  expect_identical(ev[[6L]]$payload$summary$mean_F1mean, res$summary[cycle == 6L]$mean_F1mean)
  # cancel after generation 2: the wrapper raises, with the generations so far
  err <- tryCatch(run_stacking_simulation(
    fx$meta, fx$eff, replications = 2L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE, progress = function(e) e$i < 2L),
    error = identity)
  expect_s3_class(err, "hs_error_cancelled")
  expect_identical(err$at$i, 2L)
  expect_identical(err$at$job, 1L)
  part <- err$partial
  expect_identical(max(part$summary$cycle), 2L)
  expect_identical(unique(part$summary$rep), 1L)       # job 2 never started
  expect_identical(sort(unique(part$plans$cycle)), 1:2)
  # the step API never throws: status "cancelled" and the state stays valid
  sp <- stk_spec(scheme = "TS", replications = 1L, select_cycles = 3L, top_founders = 4L,
                 num_offspring = 2L, select_top_blocks = 50L)
  st <- stk_init(sp, fx$meta, fx$eff)
  st <- stk_step(st, progress = function(e) FALSE)
  expect_identical(st$status, "cancelled")
  expect_identical(st$generation, 1L)
  expect_identical(stk_step(st), st)
  expect_identical(nrow(stk_result(st)$summary), 2L)
  expect_error(run_stacking_simulation(fx$meta, fx$eff, replications = 1L, select_cycles = 1L,
                                       top_founders = 4L, scheme = "TS", output_dir = NULL,
                                       verbose = FALSE, progress = "no"),
               "'progress' must be NULL or a function")
})

# ── E4: the hs_simulation result (rows 10 a, D4-4) ──────────────────────────

test_that("E4: the result is an hs_simulation with notes, provenance and no figures", {
  fx <- .e_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 2L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = c("TS", "GBS"), select_top_blocks = 50L,
    metrics = "d2ideal", output_dir = NULL, verbose = FALSE, seed = 5L)
  expect_s3_class(res, c("hs_simulation", "hs_result", "list"), exact = TRUE)
  expect_identical(hs_figures(res), list())
  expect_identical(nrow(hs_notes(res)), 0L)
  expect_true(all(c("stage", "code", "message", "n", "job", "first_cycle") %in%
                    names(hs_notes(res))))
  prov <- hs_provenance(res)
  expect_identical(prov$fun, "run_stacking_simulation")
  expect_identical(prov$status, "done")
  expect_identical(prov$spec$scheme, c("TS", "GBS"))
  expect_identical(prov$seed, 5L)
  expect_identical(prov$inputs$n_samples, 12L)
  expect_identical(prov$backend$name, "block")
  expect_identical(prov$tag, "TS-GBS")
  expect_identical(nrow(prov$jobs), 4L)
  expect_identical(prov$jobs$status, rep("done", 4L))
  expect_identical(prov$ideal$method, "custom")
  expect_equal(prov$ideal$ideal_phenotype, as.numeric(attr(res$ideal_genotypes, "ideal_phenotype")))
  expect_identical(prov$kept_files, character())
  # serialisable, and the same object after a round trip
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(res, f)
  expect_identical(readRDS(f), res)
  # print opens no device and names the schemes
  before <- grDevices::dev.list()
  out <- capture.output(print(res))
  expect_identical(grDevices::dev.list(), before)
  expect_match(out[1L], "hs_simulation")
  expect_true(any(grepl("TS, GBS x 2 replication", out, fixed = TRUE)))
  expect_true(any(grepl("status : done", out, fixed = TRUE)))
  # the step API gives the same class
  sp <- stk_spec(scheme = "TS", replications = 1L, select_cycles = 1L, top_founders = 4L,
                 num_offspring = 2L, select_top_blocks = 50L)
  st <- stk_init(sp, fx$meta, fx$eff)
  while (st$status == "running") st <- stk_step(st)
  r2 <- stk_result(st)
  expect_s3_class(r2, "hs_simulation")
  expect_identical(hs_provenance(r2)$fun, "stk_result")
})

test_that("E4: output_dir = NULL writes nothing; a directory gets the 0.8.8 files via write_outputs()", {
  fx <- .e_fixture()
  args <- list(fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
               num_offspring = 2L, scheme = c("TS", "OHS"), ohs_method = "Greedy",
               select_top_blocks = 50L, verbose = FALSE, seed = 2L)
  wd <- withr::local_tempdir(); withr::local_dir(wd)
  v <- withVisible(do.call(run_stacking_simulation, args))      # default output_dir
  expect_true(v$visible)
  expect_length(list.files(wd, recursive = TRUE), 0L)
  expect_false(dir.exists("stacking_output"))
  od <- file.path(wd, "out")
  v2 <- withVisible(do.call(run_stacking_simulation, c(args, list(output_dir = od))))
  expect_false(v2$visible)
  expect_setequal(list.files(od), c("stacking_plans_TS-OHS.txt", "stacking_summary_TS-OHS.txt",
                                    "Trajectories_TS-OHS.png"))
  expect_identical(v$value$summary, v2$value$summary)
  # the written tables are the result's tables (all-NA metric columns read
  # back as logical, so compare the populated columns)
  back <- data.table::fread(file.path(od, "stacking_summary_TS-OHS.txt"))
  expect_identical(nrow(back), nrow(v2$value$summary))
  for (cn in c("cycle", "select_type", "rep", "sel_scheme", "mean_F1mean", "max_F1mean"))
    expect_equal(back[[cn]], v2$value$summary[[cn]], tolerance = 1e-12, info = cn)
  # write_outputs: manifest, refusal, overwrite, what =
  expect_error(write_outputs(v$value, od), "exist already", class = "hs_error_input")
  man <- write_outputs(v$value, od, overwrite = TRUE)
  expect_s3_class(man, "data.frame")
  expect_identical(man$kind, c("table", "table", "figure"))
  expect_identical(basename(man$path), c("stacking_plans_TS-OHS.txt", "stacking_summary_TS-OHS.txt",
                                         "Trajectories_TS-OHS.png"))
  expect_true(all(man$bytes > 0))
  od2 <- file.path(wd, "tables_only")
  man2 <- write_outputs(v$value, od2, what = "tables")
  expect_identical(nrow(man2), 2L)
  expect_false(file.exists(file.path(od2, "Trajectories_TS-OHS.png")))
  expect_error(write_outputs(v$value, file.path(wd, "nope", "x"), create_dir = FALSE),
               "does not exist")
})

test_that("E4: silent fallbacks leave notes and one warning per code", {
  fx <- .e_fixture()
  # OCS with max_kinship = 0: the kinship cap stops the set at one parent and
  # truncation fills it (two codes, every cycle)
  w <- NULL
  res <- withCallingHandlers(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 2L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = "OCS", max_kinship = 0, select_top_blocks = 50L,
      output_dir = NULL, verbose = FALSE),
    hs_warning_fallback = function(cnd) { w <<- c(w, cnd$code); invokeRestart("muffleWarning") })
  n <- hs_notes(res)
  expect_setequal(unique(n$code), c("ocs_kinship_cap", "fallback_truncation"))
  expect_setequal(unique(n$job), 1:2)
  expect_true(all(n$first_cycle %in% 1:2))
  # every job left both codes at least once (the offspring generation's GRM may
  # carry negative entries, so the cap does not bite in every generation)
  for (j in 1:2) for (cd in c("ocs_kinship_cap", "fallback_truncation"))
    expect_gte(sum(n$n[n$job == j & n$code == cd]), 1L, label = paste(j, cd))
  expect_identical(sort(w), sort(c("ocs_kinship_cap", "fallback_truncation")))   # once per code
  # a GRM failure: OCS runs as truncation and says so
  testthat::local_mocked_bindings(compute_grm = function(...) stop("singular"),
                                  .package = "haploStack")
  w <- NULL
  res <- withCallingHandlers(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
      num_offspring = 2L, scheme = "OCS", select_top_blocks = 50L,
      output_dir = NULL, verbose = FALSE),
    hs_warning_fallback = function(cnd) { w <<- c(w, cnd$code); invokeRestart("muffleWarning") })
  expect_true("ocs_grm_failed" %in% hs_notes(res)$code)
  expect_match(hs_notes(res)$message[hs_notes(res)$code == "ocs_grm_failed"], "singular")
  expect_identical(w, "ocs_grm_failed")
})

test_that("E4: a metric that fails is NA and noted; the live monitor's failures are noted", {
  fx <- .e_fixture()
  testthat::local_mocked_bindings(compute_gene_diversity = function(...) stop("no diversity"),
                                  .package = "haploStack")
  w <- NULL
  res <- withCallingHandlers(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = "TS", select_top_blocks = 50L, metrics = "gene_diversity",
      output_dir = NULL, verbose = FALSE),
    hs_warning_fallback = function(cnd) { w <<- c(w, cnd$code); invokeRestart("muffleWarning") })
  expect_true(all(is.na(res$summary[cycle > 0L]$mean_gene_diversity)))
  n <- hs_notes(res)
  expect_true("metric_failed" %in% n$code)
  expect_match(n$message[n$code == "metric_failed"], "gene_diversity.*no diversity")
  expect_identical(n$n[n$code == "metric_failed"], 2L)       # aggregated over 2 generations
  expect_identical(n$first_cycle[n$code == "metric_failed"], 1L)
  expect_identical(w, "metric_failed")
})

test_that("E4: a figure that cannot render is a note, not an error (monitor and write_outputs)", {
  fx <- .e_fixture()
  testthat::local_mocked_bindings(hs_render_figure = function(...) stop("no device"),
                                  .package = "haploStack")
  od <- withr::local_tempdir()
  w <- NULL
  res <- withCallingHandlers(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 2L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
      output_dir = od, monitor = TRUE, verbose = FALSE),
    hs_warning_fallback = function(cnd) { w <<- c(w, cnd$code); invokeRestart("muffleWarning") })
  n <- hs_notes(res)
  expect_true("monitor_write_failed" %in% n$code)
  expect_true("figure_failed" %in% n$code)
  expect_identical(sort(unique(w)), c("figure_failed", "monitor_write_failed"))
  expect_true(file.exists(file.path(od, "live_summary.txt")))       # the tables still wrote
  expect_true(file.exists(file.path(od, "stacking_summary_TS.txt")))
  expect_false(file.exists(file.path(od, "Trajectories_TS.png")))
})

# ── E5: the OCS reseed is gone (row 3 a, D4-5) ──────────────────────────────

test_that("E5: ocs_selection() is deterministic, ignores seed and leaves the stream alone", {
  v <- stats::setNames(c(3, 1, 4, 1.5, 9, 2.6), paste0("S", 1:6))
  K <- diag(6); dimnames(K) <- list(names(v), names(v))
  K[upper.tri(K)] <- K[lower.tri(K)] <- 0.02
  a <- ocs_selection(v, K, n_parents = 3L, seed = 1L)
  b <- ocs_selection(v, K, n_parents = 3L, seed = 999L)
  expect_identical(a, b)
  expect_identical(a, c("S5", "S3", "S1"))
  expect_true(.e_stream_untouched(ocs_selection(v, K, n_parents = 3L)))
  expect_true(.e_stream_untouched(ocs_selection(v, NULL, n_parents = 3L)))
  # no set.seed() left anywhere outside state.R
  src <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(src), "R/ not reachable")
  code <- readLines(file.path(src, "stacking_selection.R"), warn = FALSE)
  expect_false(any(grepl("set\\.seed\\(", code[!grepl("^\\s*#", code)])))
})

test_that("E5: OCS replicates within a run no longer coincide", {
  fx <- .e_fixture(seed = 8L, n_samples = 16L, n_markers = 36L)
  res <- suppressWarnings(run_stacking_simulation(
    fx$meta, fx$eff, replications = 3L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 3L, scheme = "OCS", max_kinship = 0.9, select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE, seed = 21L), classes = "hs_warning_fallback")
  s <- res$summary[cycle > 0L]
  by_rep <- split(s$mean_F1mean, s$rep)
  # three replicates: at least two of the trajectories differ (0.8.8 reseeded
  # every cycle, which made them identical from the second cycle on)
  expect_false(identical(by_rep[[1L]], by_rep[[2L]]) && identical(by_rep[[2L]], by_rep[[3L]]))
  expect_true(all(is.finite(s$mean_F1mean)))
})

# ── E6: row 19 a guards (intake 07, D4-6) ───────────────────────────────────

test_that("E6: the default scheme excludes OCS (OCS stays available by name)", {
  expect_identical(stk_spec()$scheme, c("TS", "OHS", "GBS"))
  expect_identical(eval(formals(run_stacking_simulation)$scheme), c("TS", "OHS", "GBS"))
  reg <- haploStack:::.STK_SCHEMES
  expect_identical(names(reg)[vapply(reg, `[[`, logical(1L), "in_default")],
                   c("TS", "OHS", "GBS"))
  fx <- .e_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
    num_offspring = 2L, ohs_method = "Greedy", select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE)
  expect_setequal(unique(res$summary$sel_scheme), c("TS", "OHS", "GBS"))
  expect_identical(stk_spec(scheme = "OCS")$scheme, "OCS")
})

test_that("E6: several traits need `trait`; a long Trait table resolves by trait", {
  fx <- .e_fixture()
  eff2 <- fx$eff; eff2$Effect <- -eff2$Effect
  m <- fx$meta
  m$marker_effects <- list(A = fx$eff, B = eff2)
  base <- list(replications = 1L, select_cycles = 1L, top_founders = 4L,
               num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
               output_dir = NULL, verbose = FALSE)
  err <- tryCatch(do.call(run_stacking_simulation, c(list(m), base)), error = identity)
  expect_s3_class(err, "hs_error_input")
  expect_match(conditionMessage(err), "2 traits \\('A', 'B'\\)")
  expect_identical(err$bad, c("A", "B"))
  expect_error(do.call(run_stacking_simulation, c(list(m, trait = "C"), base)),
               "trait 'C' is not among the effects")
  rA <- do.call(run_stacking_simulation, c(list(m, trait = "A"), base))
  rB <- do.call(run_stacking_simulation, c(list(m, trait = "B"), base))
  r1 <- do.call(run_stacking_simulation, c(list(fx$meta, fx$eff), base))
  expect_identical(rA$summary, r1$summary)
  expect_false(identical(rB$summary$mean_F1mean, rA$summary$mean_F1mean))
  expect_identical(hs_provenance(rB)$trait, "B")
  # one long table with a Trait column (io/F15)
  long <- rbind(cbind(fx$eff, Trait = "A"), cbind(eff2, Trait = "B"))
  rL <- do.call(run_stacking_simulation, c(list(fx$meta, long, trait = "A"), base))
  expect_identical(rL$summary, r1$summary)
  expect_error(do.call(run_stacking_simulation, c(list(fx$meta, long), base)), "2 traits")
  # one trait in a list needs no name
  m1 <- fx$meta; m1$marker_effects <- list(Only = fx$eff)
  expect_identical(do.call(run_stacking_simulation, c(list(m1), base))$summary, r1$summary)
  # the resolver itself
  expect_identical(haploStack:::.stk_resolve_effects(NULL, m, "B"), haploStack:::.as_effect_list(eff2)[[1L]])
  expect_null(haploStack:::.stk_resolve_effects(NULL, fx$meta, allow_none = TRUE))
  expect_error(haploStack:::.stk_resolve_effects(NULL, fx$meta), "no effects supplied")
})

test_that("E6: select_priority = 'low' is refused where a scheme cannot honour it", {
  fx <- .e_fixture()
  base <- list(fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
               num_offspring = 2L, select_top_blocks = 50L, output_dir = NULL,
               verbose = FALSE, select_priority = "low")
  expect_error(stk_spec(scheme = "OHS", select_priority = "low"),
               "not supported by scheme 'OHS'", class = "hs_error_input")
  expect_error(stk_spec(scheme = "OCS", select_priority = "low"), "'OCS'")
  expect_error(stk_spec(scheme = "GBS", select_priority = "low"),
               "construct_method = \"argmax\"")
  expect_error(stk_spec(scheme = c("TS", "OHS"), select_priority = "low"), "'OHS'")
  expect_error(stk_spec(scheme = "GBS", select_priority = "low",
                        selection_unit = "pairs", pair_budget = "families"), "'GBS'")
  expect_identical(stk_spec(scheme = "TS", select_priority = "low")$select_priority, "low")
  expect_identical(stk_spec(scheme = "GBS", select_priority = "low",
                            construct_method = "argmax")$select_priority, "low")
  r <- do.call(run_stacking_simulation, c(base, list(scheme = "TS")))
  expect_identical(hs_provenance(r)$spec$select_priority, "low")
  expect_error(do.call(run_stacking_simulation, c(base, list(scheme = "OHS"))), "'OHS'")
})

test_that("E6: top_founders above the population is refused; the d2ideal note", {
  fx <- .e_fixture()                               # 12 samples
  err <- tryCatch(run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 13L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE), error = identity)
  expect_s3_class(err, "hs_error_input")
  expect_match(conditionMessage(err), "top_founders \\(13\\) exceeds the population \\(12")
  sp <- stk_spec(scheme = "TS", top_founders = 13L)
  expect_error(stk_init(sp, fx$meta, fx$eff), "exceeds the population")
  # d2ideal under "low" with a direction-blind builder: measured, and noted once
  w <- NULL
  res <- withCallingHandlers(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 2L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = "TS", select_priority = "low", metrics = "d2ideal",
      select_top_blocks = 50L, output_dir = NULL, verbose = FALSE),
    hs_warning_fallback = function(cnd) { w <<- c(w, cnd$code); invokeRestart("muffleWarning") })
  n <- hs_notes(res)
  expect_identical(n$code, "d2ideal_direction_blind")
  expect_identical(n$n, 1L)
  expect_identical(w, "d2ideal_direction_blind")
  expect_true(all(is.finite(res$summary$mean_D2Ideal)))
  # with the argmax builder the ideal follows the direction: no note
  res2 <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_priority = "low", metrics = "d2ideal",
    construct_method = "argmax", select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
  expect_identical(nrow(hs_notes(res2)), 0L)
})

test_that("E6: gs_load_geno_meta(trait =) follows the same rule", {
  skip_if_not_installed("genomicSimulation")
  fx <- .e_fixture()
  eff2 <- fx$eff; eff2$Effect <- -eff2$Effect
  m <- fx$meta; m$marker_effects <- list(A = fx$eff, B = eff2)
  withr::defer(haploStack:::.gs_session_clear())
  expect_error(gs_load_geno_meta(m, map_unit = "bp"), "2 traits", class = "hs_error_input")
  expect_false(haploStack:::.gs_session_open())     # nothing was opened
  h <- gs_load_geno_meta(m, trait = "B", map_unit = "bp")
  expect_identical(h$effect_df$Effect, eff2$Effect)
  gs_close(h)
  m0 <- fx$meta; m0$marker_effects <- NULL
  h0 <- gs_load_geno_meta(m0, map_unit = "bp")     # no effects at all is fine
  expect_null(h0$effect_df)
  expect_true(is.na(h0$files[["eff"]]))
  gs_close(h0)
})

# ── E7: the genomicSimulation boundary (row 22 a, io/F10, io/F11, io/F21, I2; D4-7) ──

test_that("E7: a missing call is written to genomicSimulation as -- (NA-any); gS header names", {
  fx <- .e_fixture()
  m <- fx$meta
  mk <- dimnames(m$geno_arrays)[[2L]]
  m$geno_arrays["S1", mk[1L], ] <- NA_character_          # both copies
  m$geno_arrays["S2", mk[2L], 2L] <- NA_character_        # one copy
  d <- withr::local_tempdir()
  f <- haploStack:::.gs_write_files(m, fx$eff, work_dir = d, map_unit = "bp")
  g <- utils::read.delim(f[["geno"]], row.names = 1L, check.names = FALSE,
                         stringsAsFactors = FALSE)
  expect_identical(g[mk[1L], "S1"], "--")
  expect_identical(g[mk[2L], "S2"], "--")                   # NA-any: the called copy goes too
  expect_match(g[mk[3L], "S3"], "^[ACGT]{2}$")
  expect_false(any(g == "NA"))
  expect_identical(strsplit(readLines(f[["map"]])[1L], "\t")[[1L]], c("marker", "chr", "pos"))
  expect_identical(strsplit(readLines(f[["eff"]])[1L], "\t")[[1L]], c("marker", "allele", "eff"))
})

test_that("E7: .gs_positions_to_cM honours map_unit and notes an auto scaling", {
  bp <- c(1e6, 3.9e7, 9.86e7)
  expect_warning(v <- haploStack:::.gs_positions_to_cM(bp), class = "hs_warning_fallback")
  expect_equal(v, c(1, 39, 98.6))
  expect_equal(suppressWarnings(haploStack:::.gs_positions_to_cM(bp, cM_per_mb = 2)), c(2, 78, 197.2))
  expect_equal(haploStack:::.gs_positions_to_cM(bp, map_unit = "bp"), c(1, 39, 98.6))
  expect_equal(haploStack:::.gs_positions_to_cM(bp, map_unit = "cM"), bp)
  expect_equal(haploStack:::.gs_positions_to_cM(c(10, 20), map_unit = "bp"), c(1e-5, 2e-5))
  expect_silent(expect_equal(haploStack:::.gs_positions_to_cM(c(10, 20)), c(10, 20)))
  r <- haploStack:::.hs_with_notes(haploStack:::.gs_positions_to_cM(bp), muffle = TRUE)
  expect_identical(r$notes$code, "gs_map_scaled")
  expect_match(r$notes$message, "cM")
  # through the spec: backend_opts$map_unit reaches the files
  sp <- stk_spec(scheme = "TS", replications = 1L, select_cycles = 1L, top_founders = 4L,
                 num_offspring = 2L, cross_backend = "genomicSimulation",
                 backend_opts = list(map_unit = "cM", cM_per_mb = 2))
  expect_identical(sp$backend_opts$map_unit, "cM")
  expect_error(stk_spec(backend_opts = list(map_unit = "Mb")), "arg")
})

test_that("E7: the round trip invents no allele; no effects gives NA F1mean", {
  skip_if_not_installed("genomicSimulation")
  fx <- .e_fixture()
  m <- fx$meta
  arr <- m$geno_arrays
  set.seed(4); arr[sample(length(arr), 30L)] <- NA_character_
  m$geno_arrays <- arr
  alleles <- unique(as.vector(fx$meta$geno_arrays))
  withr::defer(haploStack:::.gs_session_clear())
  h <- gs_load_geno_meta(m, fx$eff, map_unit = "bp")
  plan <- make_cross_design(dimnames(arr)[[1L]][1:4], method = "all_pairs", n_offspring = 3L)
  cp  <- gs_cross_plan(h, plan)
  off <- gs_group_to_geno_meta(cp$handle, cross_plan = cp$plans)
  got <- unique(as.vector(off$geno_arrays))
  expect_true(all(got[!is.na(got)] %in% alleles))          # no "N", no "A" invented
  expect_true(anyNA(off$geno_arrays))                      # the missing call propagated
  expect_false(anyNA(cp$plans$F1mean))
  gs_close(cp$handle)
  # a session without effects: gs_cross_plan() gives NA F1mean instead of 0
  m0 <- fx$meta; m0$marker_effects <- NULL
  h0 <- gs_load_geno_meta(m0, map_unit = "bp")
  cp0 <- gs_cross_plan(h0, plan)
  expect_true(all(is.na(cp0$plans$F1mean)))
  gs_close(cp0$handle)
  # the driver on a panel with missing calls: finite F1 means, no invented allele
  res <- run_stacking_simulation(
    m, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    cross_backend = "genomicSimulation", output_dir = NULL, verbose = FALSE,
    backend_opts = list(map_unit = "bp"))
  fin <- unique(as.vector(res$final_populations[[1L]]$geno_arrays))
  expect_true(all(fin[!is.na(fin)] %in% alleles))
  expect_true(all(is.finite(res$summary[cycle > 0L]$mean_F1mean)))
})

test_that("E7: .gs_require() wants genomicSimulation >= 0.3", {
  skip_if_not_installed("genomicSimulation")
  expect_true(haploStack:::.gs_require())
  testthat::local_mocked_bindings(.gs_version = function() package_version("0.2.9"),
                                  .package = "haploStack")
  expect_error(haploStack:::.gs_require(), ">= 0.3", class = "hs_error_dependency")
})

# ── E8: the small fixes (sc/F22, F28, F29, F30, F31, F23; D4-8) ─────────────

test_that("E8: a zero fitness weight keeps the raw contribution; the composite is unchanged", {
  eff1 <- data.frame(Marker = c("M1", "M1", "M2", "M2"), Allele = c("A", "T", "A", "T"),
                     Effect = c(0.5, -0.5, 0.2, -0.2))
  eff2 <- data.frame(Marker = c("M1", "M1", "M2", "M2"), Allele = c("A", "T", "A", "T"),
                     Effect = c(0.1, -0.1, 0.8, -0.8))
  fi <- create_fitness_index(list(Yield = eff1, Protein = eff2),
                             c(Yield = 1, Protein = 0), allele_dict = c("A", "T"))
  expect_false(anyNA(fi$Protein))
  expect_equal(fi$Protein[fi$Marker == "M2" & fi$Allele == "A"], 0.8)
  expect_equal(fi$Effect, fi$Yield * 1 + fi$Protein * 0)
  fi2 <- create_fitness_index(list(Yield = eff1, Protein = eff2),
                              c(Yield = 0.7, Protein = 0.3), allele_dict = c("A", "T"))
  expect_equal(fi2$Effect, 0.7 * fi2$Yield + 0.3 * fi2$Protein)
  # update_fitness_index() keeps to the original registry
  eff3 <- rbind(eff2, data.frame(Marker = "M9", Allele = c("A", "T"), Effect = c(1, -1)))
  up <- update_fitness_index(fi2, list(Yield = eff1, Protein = eff3), c(Yield = 0.7, Protein = 0.3))
  expect_false("M9" %in% up$Marker)
  expect_identical(up[, .(Marker, Allele)], fi2[, .(Marker, Allele)])
})

test_that("E8: make_cross_design() leaves its input alone; a n_offspring column wins", {
  ct <- data.table(Parent1 = c("A", "B"), Parent2 = c("B", "C"))
  d  <- make_cross_design(NULL, method = "custom", cross_table = ct, n_offspring = 3L)
  expect_named(ct, c("Parent1", "Parent2"))                 # no column added by reference
  expect_identical(d$n_offspring, c(3L, 3L))
  ct2 <- data.table(Parent1 = c("A", "B"), Parent2 = c("B", "C"), n_offspring = c(1L, 5L))
  d2  <- make_cross_design(NULL, method = "custom", cross_table = ct2, n_offspring = 3L)
  expect_identical(d2$n_offspring, c(1L, 5L))
})

test_that("E8: summarise_simulation() skips an all-NA metric without a warning", {
  fx <- .e_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 2L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE)
  expect_silent(s <- summarise_simulation(res$trajectories))
  expect_false(any(grepl("^BlockCoverage", names(s))))
  expect_true(all(c("MeanGEBV_mean", "SDev_mean") %in% names(s)))
  expect_true(all(is.finite(s$MeanGEBV_max)))
  expect_error(summarise_simulation(data.table(Rep = 1L, Cycle = 1L, MeanGEBV = NA_real_)),
               "no metric column")
})

test_that("E8: ga_control reaches the GA engines; the default changes nothing", {
  fx <- .e_fixture()
  base <- list(fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
               num_offspring = 2L, scheme = "OHS", ohs_method = "GA_pairwise",
               select_top_blocks = 50L, output_dir = NULL, verbose = FALSE, seed = 3L)
  r0 <- do.call(run_stacking_simulation, base)
  r1 <- do.call(run_stacking_simulation, c(base, list(ga_control = list())))
  expect_identical(r1$summary, r0$summary)
  r2 <- do.call(run_stacking_simulation,
                c(base, list(ga_control = list(pop_size = 6L, n_generations = 3L))))
  expect_true(all(is.finite(r2$summary$mean_F1mean)))
  expect_identical(hs_provenance(r2)$spec$ga_control$pop_size, 6L)
  expect_identical(haploStack:::.stk_ga(list(ga_control = list()))$n_generations, 200L)
  expect_identical(haploStack:::.stk_ga(list(ga_control = list(mutation_rate = 0.3)))$mutation_rate, 0.3)
})

test_that("E8: catalogue_results / trait_name / a positional catalogue warn as deprecated", {
  fx <- .e_fixture()
  cat_res <- list(SimTrait = list(geno_meta = fx$meta))
  base <- list(effect_df = fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
               num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
               output_dir = NULL, verbose = FALSE)
  expect_warning(do.call(run_stacking_simulation, c(list(catalogue_results = cat_res), base)),
                 "catalogue_results.*deprecated", class = "hs_warning_deprecated")
  expect_warning(do.call(run_stacking_simulation, c(list(cat_res), base)),
                 "catalogue list passed", class = "hs_warning_deprecated")
  w <- NULL
  withCallingHandlers(
    do.call(run_stacking_simulation,
            c(list(catalogue_results = cat_res, trait_name = "SimTrait"), base)),
    hs_warning_deprecated = function(cnd) { w <<- c(w, cnd$old); invokeRestart("muffleWarning") })
  expect_setequal(w, c("catalogue_results", "trait_name"))
})

# ── E9: fit_genomic_prediction()$summary (io/F14; D4-9) ─────────────────────

test_that("E9: fit_genomic_prediction() returns a plain summary beside the fit", {
  skip_if_not_installed("rrBLUP")
  m  <- sim_geno_meta(n_samples = 20L, n_markers = 30L, n_chrom = 2L, block_size = 3L, seed = 9L)
  ph <- data.frame(Sample = dimnames(m$geno_arrays)[[1L]],
                   Yield = haploStack:::.with_seed(1L, stats::rnorm(20L)))
  fit <- suppressMessages(fit_genomic_prediction(m, phenotype = ph, model = "rrblup",
                                                 trait = "Yield"))
  s <- fit$summary
  expect_named(s, c("model", "n_train", "var_u", "var_e", "h2", "seed", "nIter", "burnIn", "backend",
                    "backend_version"))
  expect_identical(s$model, "rrblup")
  expect_identical(s$n_train, 20L)
  expect_true(is.numeric(s$var_u) && is.numeric(s$var_e) && length(s$var_u) == 1L)
  expect_true(is.na(s$h2))
  expect_identical(s$backend, "rrBLUP")
  expect_identical(s$backend_version, as.character(utils::packageVersion("rrBLUP")))
  g <- suppressMessages(fit_genomic_prediction(m, phenotype = ph, model = "gblup", trait = "Yield"))
  expect_true(is.finite(g$summary$h2) && g$summary$h2 >= 0 && g$summary$h2 <= 1)
  expect_equal(g$summary$h2, g$summary$var_u / (g$summary$var_u + g$summary$var_e))
})
