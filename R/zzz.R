# =============================================================================
# zzz.R — package load hooks + data.table hardening
# =============================================================================

.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")
  packageStartupMessage(pkgname, " ", version,
                        "  —  haplotype block analysis & stacking simulation")
}

# ---------------------------------------------------------------------------
# Package options
#
# haploStack.use_rcpp — global switch for ALL compiled kernels (src/*.cpp):
# the local-GEBV kernel (src/local_gebv.cpp) and the marker allele-count / dosage
# kernels (src/marker_stats.cpp). Default TRUE: each call site (gated via
# .hs_use_rcpp() in R/catalogue_engine.R) uses its Rcpp kernel when loadable and
# falls back to an identical pure-R path otherwise. Set
# options(haploStack.use_rcpp = FALSE) to force the pure-R paths (e.g. for
# debugging or A/B benchmarking). Only set the default when the user has not
# already configured it.
#
# haploStack.fill_na — the placeholder written in place of an unobserved allele
# when a block haplotype signature is assembled (.compute_block_signatures(),
# write_block_allele_sheet(), count_block_subtypes()). Default "-", the
# canonical missing-allele character of DEFAULT_ENCODING_DICT. Must not be
# empty: a signature is a positional string, so a zero-width placeholder
# shortens it and lets distinct
# allele combinations collide — over (M1, M2, M3), (A, T, -) and
# (A, -, T) would both paste to "AT". Read via .hs_fill_na() (R/alphabet.R);
# every call site also takes a per-call `fill_na` argument. Only set the default
# when the user has not already configured it.
#
# haploStack.sample_order — the broadcast row order read by every figure with a
# sample axis when its own `sample_order` argument is NULL (plot_order.R). Holds
# an hs_sample_clustering (from cluster_samples()) or a character vector of
# sample IDs. Deliberately UNSET by default — there is no sensible package-wide
# row order, and the fallback is each figure's own default (alphabetical for
# plot_haplotype_alignment(), family-history order for
# plot_localgebv_alignment()). Set it with use_sample_order(), not options():
# the setter validates, so a bad order fails at the call the user can see rather
# than inside the next figure. Precedence everywhere: per-call argument > this
# option > the figure's default, matching haploStack.fill_na.
#
# haploStack.missing_chars — extra missing-allele characters, registered per
# source format for the session (list(nucleotide = "#", ...)). Read via
# .extra_missing_chars() inside .missing_tokens(), so one registration reaches
# every consumer at once. Deliberately UNSET by default and deliberately not
# persisted: it is a statement about a file being read, not about the package,
# and a character that means "absent" in one vocabulary can be a genotype in
# another. Set it with use_missing_chars(), not options(): the setter validates
# (single character, never a real allele of that mode), so a bad registration
# fails at the call the user can see rather than inside the next read. A
# container normalised while one was active records the vocabulary in
# metadata$missing_tokens.
# ---------------------------------------------------------------------------
.onLoad <- function(libname, pkgname) {
  if (is.null(getOption("haploStack.use_rcpp")))
    options(haploStack.use_rcpp = TRUE)
  if (is.null(getOption("haploStack.fill_na")))
    options(haploStack.fill_na = HS_FILL_NA_DEFAULT)
  invisible()
}

# ---------------------------------------------------------------------------
# data.table awareness
#
# `import(data.table)` (see haploStack-package.R) already marks this namespace
# data.table-aware via cedta(); setting `.datatable.aware` makes it explicit and
# unconditional, including for `:=`/`.SD` calls evaluated inside future /
# future.apply workers in run_stacking_simulation().
# ---------------------------------------------------------------------------
.datatable.aware <- TRUE

# ---------------------------------------------------------------------------
# Global-variable registration
#
# data.table's special read-only symbols and the column names referenced bare in
# non-standard evaluation (`dt[, col := ...]`, `by = .(col)`) are not lexical
# bindings, so R CMD check reports "no visible binding for global variable".
# Declaring them here silences those NOTEs (which --as-cran / strict CI promote
# to failures). The list is the union of bare symbols across R/.
# ---------------------------------------------------------------------------
utils::globalVariables(c(
  # data.table special symbols
  ".", ".N", ".I", ".SD", ".GRP", ".BY", ".EACHI", ".NGRP",
  # bare column / aggregation names used in NSE
  "Alt", "Dose", "IsAlt", "K", "N", "Ref",
  "Allele", "AlleleCount", "AlleleFreq", "Block", "Block_chr", "Block_end",
  "Block_id", "Block_length", "Block_start", "EffKey", "Num_SNPs",
  "Num_SNPs_in_Block",
  "CFreq", "Child", "Chrom", "Chrom_", "Count", "CrossID", "Cross_ID", "Cycle",
  "Effect", "End", "F1count", "F1mean", "FamID", "Family", "Freq", "GEBV",
  "GenoEffect", "GenomicPos", "NewID", "OldID",
  "HapEffect", "HapN", "Haplotype", "Haplotype_idx", "IBD", "Label", "MAF",
  "Marker", "MarkerTotal", "Method", "P1_Allele", "P2_Allele", "Parent",
  "Parent1", "Parent2", "ParentHap", "Pos", "Pos_", "Progeny", "ProgenyHap",
  "Sample", "SampleHap", "SelFreq", "Similarity", "Start", "Subtype", "TP",
  "TargetFreq",
  "Trait", "Value", "accuracy", "allele", "block", "block_id", "candidate",
  "chr_max", "chromosome", "code", "comparable", "contrib", "f1", "j",
  "marker", "marker_same", "mcc",
  "mid", "n_offspring", "precision", "recall", "sample", "score_val",
  "signature", "sim",
  "specificity", "tkey", "x_center", "x_pos", "xmin", "y_label",
  # signal-marker selection (catalogue_similarity.R / catalogue_viz.R)
  "N_Signal", "Signal_Marker", "Signal_Score", "block_highlighted",
  "block_n_markers", "block_score", "highlighted", "is_signal", "score_column",
  "signal_rank", "tie_key", "xcenter", "xmax",
  # subtype search / carriage counting (block_subtypes.R)
  "Carried", "Copies", "Frac", "Target", "Total", "nCarried", "nTarget",
  # favourable-block counting by local-GEBV sign (stacking_evaluate.R)
  "FID", "GID", "IID", "Preferred", "Sign", "gen_i", "is_off",
  "nBlock", "nNegative", "nPositive", "nPreferred",
  # block-subtype selection value columns (catalogue_analysis.R: the raw effect
  # S and the selection base Vbase; reused by .resolve_top_subtypes())
  "S", "Vbase",
  # data.table update-join RHS symbols (`x[i, on=, col := i.col]`). These are
  # supplied by the join, never lexically bound, so every use is a check NOTE.
  # Only the four used across R/ are listed; add here when a new one appears.
  "i.Block", "i.Chrom", "i.Effect", "i.Pos",
  "i.FID", "i.GID", "i.Label"
))
