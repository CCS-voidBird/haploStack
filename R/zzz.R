# =============================================================================
# zzz.R — package load hooks + data.table hardening
# =============================================================================

# The second line is how a coding agent meeting the package in a clean
# environment finds the usage guide shipped in inst/skills (R/agent_guide.R).
.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")
  packageStartupMessage(pkgname, " ", version,
                        "  \u2014  haplotype block analysis & stacking simulation\n",
                        "Usage guide (people and AI agents): hs_agent_guide()")
}

# ---------------------------------------------------------------------------
# Package options
#
# Six settings, one rule (0.9.0, design framework F1): every exported function
# whose result depends on a setting takes it as an argument (missing_chars=,
# sample_order=, fill_na=, strict_missing=, st_formats=; use_rcpp stays a
# process option),
# and the option is only the default when the argument is NULL. They are
# resolved once, at the entry point, through .hs_setting() (R/state.R) — the
# registry .HS_SETTINGS there lists option name, default and validator — and
# passed down as plain arguments. No other file reads getOption("haploStack.")
# (test-090-hygiene.R). A value set with options() by hand is validated on
# read: a wrong one is ignored where 0.8.8 ignored it, never obeyed silently.
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
# .extra_missing_chars() ONLY where a container is built or repaired (the
# readers, build_geno_meta(), new_geno_meta(), repair_geno_meta(), each of which
# also takes `missing_chars=`); the container then records the vocabulary it
# was normalised against in metadata$missing_tokens, and that stamp is what
# every consumer goes by (.meta_extra_missing()). Deliberately UNSET by default
# and deliberately not persisted: it is a statement about a file being read,
# not about the package, and a character that means "absent" in one vocabulary
# can be a genotype in another. Set it with use_missing_chars(), not options():
# the setter validates (single character, never a real allele of that mode), so
# a bad registration fails at the call the user can see rather than inside the
# next read.
#
# haploStack.strict_missing — FALSE by default. TRUE turns the "container still
# carries a missing token" warning of is_geno_meta() / load_geno_meta() into an
# error (.notify_missing_tokens(), R/geno_meta_missing.R); both functions take
# `strict_missing=` per call. Never set by the package.
#
# haploStack.st_formats — FALSE by default. SelectionTools is being retired;
# everything that reads or writes one of its own formats (st_convert(),
# st_ld_blocks(), st_ld_blocks_files(), R/st_formats.R) refuses to run until
# this is TRUE — use_st_formats(TRUE) for the session, or `st_formats = TRUE`
# on the call. Never set by the package.
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
  "Parent1", "Parent2", "ParentHap", "Pos", "Pos_", "Pos2", "Progeny",
  "ProgenyHap",
  "Sample", "SampleHap", "SelFreq", "Similarity", "Start", "Subtype", "TP",
  "TargetFreq",
  "Trait", "Value", "accuracy", "allele", "block", "block_id", "candidate",
  "chr_max", "chromosome", "code", "comparable", "contrib", "f1", "j",
  "marker", "marker_same", "mcc",
  "mid", "n_offspring", "precision", "recall", "sample", "score_val",
  "signature", "sim", "srow", "mrow",
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
  # block annotation (annotation_read.R / annotation_compute.R, thread #033)
  "gidx", "bidx", "midx", "ridx", "qidx", "MChrom", "qs", "qe", "ms", "me",
  "chr_key", "chrom", "start", "end", "name", "strand", "val",
  "gene_id", "gene_name", "biotype", "gene_count", "gene_ids", "gene_names",
  "qtl_count", "qtl_names", "feature_type", "feature_id", "feature_name",
  "feature_biotype", "feature_start", "feature_end", "feature_strand",
  "gene_symbol", "gene_description", "go_terms", "interpro_domains",
  "gene_start", "gene_end", "gene_strand", "ensembl_browser_url",
  "gene_stable_id", "description", "xref_id",
  # circos meta plot / overlay (catalogue_circos_meta.R,
  # catalogue_annotation_overlay.R, thread #033 stage 6)
  ".val", ".plotval", "Block_Var", "Start_Pos", "End_Pos", "hi", "Name",
  # data.table update-join RHS symbols (`x[i, on=, col := i.col]`). These are
  # supplied by the join, never lexically bound, so every use is a check NOTE.
  # Only the ones used across R/ are listed; add here when a new one appears.
  "i.Block", "i.Chrom", "i.Effect", "i.Pos",
  "i.FID", "i.GID", "i.Label",
  "i.gene_count", "i.gene_ids", "i.gene_names", "i.qtl_count", "i.qtl_names"
))
