# =============================================================================
# genomic_prediction.R
#
# Bridge between the GenoMeta container and the standard genomic-prediction
# engines (rrBLUP, BGLR). Estimates per-marker allele effects and per-sample
# GEBVs from phenotypes and feeds the results back into a GenoMeta in the exact
# shapes the rest of the package consumes:
#
#   * marker_effects — a named list with one (Marker, Allele, Effect) table per
#     trait (the list name is the trait); drop-in for the effect_list of
#     compute_profile_localgebv() / run_haplotype_catalogue().
#   * a GEBV "phenotype" table (Sample, <gebv col>) — drop-in for
#     attach_phenotype().
#
# The genotype lives in geno_meta$geno_arrays as a 3-D *haplo-wise* character
# array [samples x markers x ploidy]; rrBLUP/BGLR need a numeric per-marker
# dosage matrix [individuals x markers]. The transform reuses geno_meta_views()
# so all encodings are normalised through one code path.
#
#   .gp_allele_dosage()       — internal: haplo-wise array -> {0,1,2} dosage M
#   .gp_resolve_pheno()       — internal: resolve + align a trait vector
#   .gp_gblup_marker_effects()— internal: back-solve marker effects from GBLUP
#   fit_genomic_prediction()  — exported: fit rrblup / gblup / bayesc
#   predict_gebv()            — exported: marker effects -> whole-genome GEBV,
#                               appended to the GenoMeta phenotype for plotting
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: haplo-wise allele array -> numeric dosage matrix
# ---------------------------------------------------------------------------

#' Build a numeric alt-allele dosage matrix from a GenoMeta
#'
#' Converts the haplo-wise genotype stored in \code{geno_meta} into the numeric
#' \code{[samples x markers]} dosage matrix required by rrBLUP / BGLR. Each cell
#' counts the copies of the marker's \emph{alternate} (minor) allele carried by a
#' sample, giving integers in \code{{0, 1, 2}} for a diploid.
#'
#' Built directly from \code{geno_meta$geno_arrays} so the matrix columns and the
#' returned \code{allele_map} carry the \strong{canonical marker names} from the
#' map (the \code{dimnames} of the array). It deliberately does \emph{not} route
#' through \code{\link{geno_meta_views}}, whose wide allele table passes marker
#' names through \code{as.data.frame()} and would mangle non-syntactic IDs
#' (e.g. \code{"1H-123"} -> \code{"X1H.123"}).
#'
#' Markers with more than two observed alleles keep their two most frequent
#' alleles (the rest are treated as reference) and trigger a \code{warning()}.
#' Missing dosages are mean-imputed per marker when \code{impute = "mean"}.
#'
#' @param geno_meta A valid \code{GenoMeta} list with a diploid allele encoding
#'   (\code{"selectiontools"} or \code{"nucleotide"}).
#' @param impute    \code{"mean"} (default) mean-imputes missing dosages per
#'   marker; \co\code{"none"} leaves them as \code{NA}.
#' @return A list with elements \code{M} (numeric matrix, rows = samples and
#'   cols = markers in \code{geno_arrays} order, names preserved verbatim),
#'   \code{allele_map} (data.frame: \code{Marker}, \code{Ref}, \code{Alt}) and
#'   \code{samples}.
#' @keywords internal
.gp_allele_dosage <- function(geno_meta, impute = c("mean", "none")) {
  impute <- match.arg(impute)

  enc <- geno_meta$metadata$encoding
  if (is.null(enc))
    enc <- detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))
  if (!enc %in% c("selectiontools", "nucleotide"))
    stop(".gp_allele_dosage: marker-effect estimation requires a diploid ",
         "allele encoding ('selectiontools' or 'nucleotide'); got '",
         enc, "'. Recode with recode_geno_meta() first.")

  # Operate on the 3-D array directly: dimnames carry the canonical marker names
  # from the map, so column/allele-map names match geno_meta$map_info$Marker.
  arr     <- geno_meta$geno_arrays
  samples <- dimnames(arr)[[1L]]
  markers <- dimnames(arr)[[2L]]
  n_samp  <- length(samples)
  n_mark  <- length(markers)
  ploidy  <- dim(arr)[3L]

  # Ref/alt detection (blanks stripped, ties A->Z) via the allele-count kernel:
  # ref = most frequent observed allele, alt = 2nd most frequent. Detection is in
  # R; the kernel only counts codes.
  #
  # `na_vals` is defence-in-depth (FR7), not the mechanism. A container built by
  # 0.6.0 arrives normalised, so the token set is empty of anything present and
  # this is a no-op. It exists for the container that did not come through the
  # portal — a pre-0.6.0 .rds, or a hand-assembled list — where without it the
  # missing code ranks as an allele and the marker's effect is estimated for
  # missingness itself (INDEX discrepancy 23).
  miss   <- .missing_tokens(enc)
  counts <- .marker_allele_counts(arr, na_vals = miss, strip_blank = TRUE)
  levs    <- colnames(counts)
  ref     <- character(n_mark)
  alt     <- character(n_mark)
  n_multi <- 0L
  for (j in seq_len(n_mark)) {
    cj  <- counts[j, ]
    obs <- which(cj > 0L)
    if (!length(obs)) { alt[j] <- NA_character_; next } # all missing at marker
    ord    <- obs[order(-cj[obs])]                      # count desc; ties A->Z
    ref[j] <- levs[ord[1L]]
    alt[j] <- if (length(ord) >= 2L) levs[ord[2L]] else NA_character_
    if (length(ord) > 2L) n_multi <- n_multi + 1L
  }

  # Alt-allele dosage; NA only when a cell is FULLY missing (na_if_any = FALSE),
  # blanks treated as observed non-matches (strip_blank = FALSE) — both preserving
  # the original per-marker semantics. Monomorphic markers (alt NA) -> 0 copies.
  # `na_vals` is the same defence-in-depth as above: a no-op on a normalised
  # container, and on a legacy one it stops the missing code being counted as an
  # observed non-match of the alt allele.
  M <- .dosage_matrix(arr, alt, na_if_any = FALSE, na_vals = miss,
                      strip_blank = FALSE)
  M <- matrix(as.numeric(M), n_samp, n_mark, dimnames = list(samples, markers))
  # Markers with no non-blank observed allele stay all-NA (matches the original
  # `next`, which left the NA_real_-initialised column untouched).
  empty_mark <- rowSums(counts) == 0L
  if (any(empty_mark)) M[, empty_mark] <- NA_real_

  if (n_multi > 0L)
    warning(".gp_allele_dosage: ", n_multi, " marker(s) have >2 alleles; ",
            "keeping the two most frequent and treating the rest as reference.")

  if (impute == "mean" && anyNA(M)) {
    for (j in seq_len(n_mark)) {
      na_j <- is.na(M[, j])
      if (any(na_j)) {
        mu <- mean(M[!na_j, j])
        M[na_j, j] <- if (is.nan(mu)) 0 else mu
      }
    }
  }

  allele_map <- data.frame(Marker = markers, Ref = ref, Alt = alt,
                           stringsAsFactors = FALSE)
  list(M = M, allele_map = allele_map, samples = samples)
}


# ---------------------------------------------------------------------------
# Internal: resolve and align a phenotype/trait vector
# ---------------------------------------------------------------------------

#' Resolve a trait vector aligned to a set of samples
#'
#' Picks the phenotype source (explicit \code{phenotype} overrides
#' \code{geno_meta$phenotype}), selects the trait column and returns a numeric
#' vector aligned by name to \code{sample_ids} (\code{NA} where a sample has no
#' record — these become the held-out / non-training rows).
#'
#' @param geno_meta  A valid \code{GenoMeta} list.
#' @param phenotype  Optional phenotype data.frame (\code{Sample} + trait
#'   columns) or file path; \code{NULL} uses \code{geno_meta$phenotype}.
#' @param trait      Trait column name; \code{NULL} picks the first numeric
#'   non-\code{Sample} column.
#' @param sample_ids Character vector of sample IDs to align to.
#' @return A list with \code{y} (named numeric vector) and \code{trait} (the
#'   resolved trait name).
#' @keywords internal
.gp_resolve_pheno <- function(geno_meta, phenotype, trait, sample_ids) {
  if (is.null(phenotype)) phenotype <- geno_meta$phenotype
  if (is.null(phenotype))
    stop(".gp_resolve_pheno: no phenotype supplied — pass 'phenotype' or ",
         "attach one to the GenoMeta with attach_phenotype().")
  if (is.character(phenotype) && length(phenotype) == 1L)
    phenotype <- read_phenotype(phenotype)
  if (!is.data.frame(phenotype))
    stop(".gp_resolve_pheno: 'phenotype' must be a data.frame or a file path.")
  phenotype <- as.data.frame(phenotype, stringsAsFactors = FALSE)

  sample_col <- if ("Sample" %in% names(phenotype)) "Sample" else names(phenotype)[1L]
  ids        <- as.character(phenotype[[sample_col]])
  trait_cols <- setdiff(names(phenotype), sample_col)

  if (is.null(trait)) {
    is_num <- vapply(phenotype[trait_cols], function(x)
      is.numeric(x) || all(!is.na(suppressWarnings(as.numeric(x)))), logical(1L))
    if (!any(is_num))
      stop(".gp_resolve_pheno: no numeric trait column found in the phenotype.")
    trait <- trait_cols[is_num][1L]
  } else if (!trait %in% trait_cols) {
    stop(".gp_resolve_pheno: trait '", trait, "' not found. Available: ",
         paste(trait_cols, collapse = ", "), ".")
  }

  vals <- suppressWarnings(as.numeric(phenotype[[trait]]))
  y    <- stats::setNames(vals[match(sample_ids, ids)], sample_ids)
  if (all(is.na(y)))
    stop(".gp_resolve_pheno: trait '", trait, "' has no values for any ",
         "genotype sample.")
  list(y = y, trait = trait)
}


# ---------------------------------------------------------------------------
# Internal: back-solve marker effects from a GBLUP solution
# ---------------------------------------------------------------------------

#' Back-solve per-marker effects from a GBLUP genomic value vector
#'
#' Recovers per-marker effects from the GBLUP solution using the RR-BLUP /
#' SNP-BLUP equivalence \code{u = t(Z) \%*\% solve(Z \%*\% t(Z)) \%*\% g}, where
#' \code{Z} is the column-centred marker matrix and \code{g} the per-sample
#' genomic values. A small ridge is added to \code{Z \%*\% t(Z)} for numerical
#' stability. The result is the minimum-norm marker-effect vector consistent with
#' \code{g = Z \%*\% u}; it is an approximation (GBLUP does not estimate marker
#' effects directly).
#'
#' @param Z A column-centred numeric marker matrix \code{[samples x markers]}.
#' @param g Named numeric vector of genomic values, aligned to \code{rownames(Z)}.
#' @return A numeric vector of marker effects, one per column of \code{Z}.
#' @keywords internal
.gp_gblup_marker_effects <- function(Z, g) {
  g   <- as.numeric(g[rownames(Z)])
  ZZt <- tcrossprod(Z)
  ridge <- 1e-6 * mean(diag(ZZt))
  diag(ZZt) <- diag(ZZt) + ridge
  as.numeric(crossprod(Z, solve(ZZt, g)))
}


# ---------------------------------------------------------------------------
# Exported: fit genomic-prediction model
# ---------------------------------------------------------------------------

#' Estimate marker effects and GEBVs from a GenoMeta with rrBLUP / BGLR
#'
#' Connects the \code{GenoMeta} container to the standard genomic-prediction
#' engines. It transforms the haplo-wise genotype array into a numeric dosage
#' matrix, fits the requested model against a phenotype, and returns the
#' estimated per-marker allele effects and per-sample GEBVs \strong{already
#' shaped} for re-use elsewhere in the package — as \code{marker_effects} (a
#' named list of per-trait tables, the effect_list shape consumed by
#' \code{\link{compute_profile_localgebv}} and
#' \code{\link{run_haplotype_catalogue}}) and as a GEBV "phenotype" table (for
#' \code{\link{attach_phenotype}}). A copy of the input \code{GenoMeta} with both
#' attached and revalidated is returned in \code{$geno_meta}; the fitted effects
#' are merged into any effects already stored there (one table per trait).
#'
#' Three models are available:
#' \describe{
#'   \item{\code{"rrblup"}}{Ridge-regression BLUP (RR-BLUP) via
#'     \code{rrBLUP::mixed.solve()}; estimates marker effects directly.}
#'   \item{\code{"gblup"}}{Genomic BLUP via \code{rrBLUP::A.mat()} +
#'     \code{rrBLUP::kin.blup()}; per-sample GEBVs are exact, and marker effects
#'     are back-solved from the GBLUP solution (see
#'     \code{.gp_gblup_marker_effects}).}
#'   \item{\code{"bayesc"}}{BayesC via \code{BGLR::BGLR()}.}
#' }
#'
#' rrBLUP (for \code{"rrblup"}/\code{"gblup"}) and BGLR (for \code{"bayesc"}) are
#' optional dependencies; install the relevant one before calling. The genotype
#' must use a diploid allele encoding (\code{"selectiontools"} or
#' \code{"nucleotide"}); recode with \code{\link{recode_geno_meta}} if needed.
#'
#' @param geno_meta    A valid \code{GenoMeta} list.
#' @param phenotype    Optional phenotype data.frame (\code{Sample} + trait
#'   columns) or file path. When \code{NULL} (default),
#'   \code{geno_meta$phenotype} is used. Samples without a phenotype value are
#'   held out of training and still receive a predicted GEBV. An external
#'   phenotype is also attached to the returned \code{$geno_meta} (restricted to
#'   the individuals shared with the genotype), so the observed trait values sit
#'   alongside the GEBV column.
#' @param model        One of \code{"rrblup"} (default), \code{"gblup"},
#'   \code{"bayesc"}.
#' @param trait        Trait column to model. \code{NULL} (default) uses the
#'   first numeric trait column.
#' @param gebv_colname Name for the GEBV column in the returned phenotype table.
#'   \code{NULL} (default) uses \code{paste0("GEBV_", trait)}.
#' @param impute       Missing-dosage handling passed to the dosage transform:
#'   \code{"mean"} (default) or \co\code{"none"}.
#' @param nIter,burnIn BGLR sampler length and burn-in (\code{"bayesc"} only).
#' @param verbose      Logical; print backend progress. Default \code{FALSE}.
#' @param ...          Extra arguments forwarded to the backend
#'   (\code{mixed.solve} / \code{kin.blup} / \code{BGLR}).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{marker_effects}}{Named list with one element — the fitted
#'     \code{trait} — holding a data.frame \code{Marker, Allele, Effect} with two
#'     rows per polymorphic marker (alternate allele = estimated effect,
#'     reference allele = 0). This is the effect_list shape; ready to pass to
#'     \code{\link{run_haplotype_catalogue}} or assign to
#'     \code{geno_meta$marker_effects}.}
#'   \item{\code{gebv}}{data.frame \code{Sample, <gebv_colname>}. Ready for
#'     \code{\link{attach_phenotype}}.}
#'   \item{\code{model}}{The fitted backend object.}
#'   \item{\code{geno_meta}}{A copy of the input with \code{marker_effects} and a
#'     phenotype holding the modelled trait values (an external phenotype is
#'     attached for the shared individuals) plus the GEBV column, revalidated
#'     with \code{is_geno_meta()}.}
#'   \item{\code{allele_map}}{data.frame \code{Marker, Ref, Alt} used for coding.}
#'   \item{\code{trait}}{The resolved trait name.}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- sim_geno_meta(40, 60, block_size = 5L, seed = 1L)
#' ph   <- data.frame(Sample = dimnames(meta$geno_arrays)[[1]],
#'                    Yield  = rnorm(40))
#' res  <- fit_genomic_prediction(meta, ph, model = "rrblup")
#' head(res$marker_effects$Yield)             # one table per trait
#' head(res$gebv)
#' compute_profile_localgebv(res$geno_meta)   # uses the fitted effects
#' run_haplotype_catalogue(res$geno_meta)     # effects detected from GenoMeta
#' }
#'
#' @seealso \code{\link{compute_profile_localgebv}}, \code{\link{attach_phenotype}},
#'   \code{\link{recode_geno_meta}}, \code{\link{geno_meta_views}}
#' @export
fit_genomic_prediction <- function(geno_meta,
                                   phenotype    = NULL,
                                   model        = c("rrblup", "gblup", "bayesc"),
                                   trait        = NULL,
                                   gebv_colname = NULL,
                                   impute       = c("mean", "none"),
                                   nIter        = 12000L,
                                   burnIn       = 2000L,
                                   verbose      = FALSE,
                                   ...) {
  model  <- match.arg(model)
  impute <- match.arg(impute)
  if (!is_geno_meta(geno_meta))
    stop("fit_genomic_prediction: 'geno_meta' is not a valid GenoMeta.")

  backend <- if (model == "bayesc") "BGLR" else "rrBLUP"
  if (!requireNamespace(backend, quietly = TRUE))
    stop("fit_genomic_prediction: model '", model, "' requires the '", backend,
         "' package. Install it with install.packages(\"", backend, "\").")

  # Genotype -> dosage matrix; phenotype -> aligned trait vector.
  d   <- .gp_allele_dosage(geno_meta, impute = impute)
  M   <- d$M
  amap <- d$allele_map
  rp  <- .gp_resolve_pheno(geno_meta, phenotype, trait, d$samples)
  y   <- rp$y
  trait <- rp$trait
  if (is.null(gebv_colname)) gebv_colname <- paste0("GEBV_", trait)

  train <- !is.na(y)
  if (sum(train) < 2L)
    stop("fit_genomic_prediction: need >= 2 phenotyped samples to train; got ",
         sum(train), ".")

  # ---- Fit ----------------------------------------------------------------
  if (model == "rrblup") {
    fit  <- rrBLUP::mixed.solve(y = y[train], Z = M[train, , drop = FALSE], ...)
    u    <- as.numeric(fit$u)
    gebv <- as.numeric(M %*% u)

  } else if (model == "gblup") {
    Mc <- M - 1                                   # {-1,0,1} for A.mat()
    A  <- rrBLUP::A.mat(Mc)
    dat <- data.frame(Sample = d$samples, y = unname(y), stringsAsFactors = FALSE)
    fit <- rrBLUP::kin.blup(data = dat, geno = "Sample", pheno = "y", K = A, ...)
    g    <- fit$g[d$samples]
    gebv <- as.numeric(g)
    Zc   <- scale(Mc, center = TRUE, scale = FALSE)
    u    <- .gp_gblup_marker_effects(Zc, stats::setNames(gebv, d$samples))

  } else {                                        # bayesc
    eta <- list(list(X = M[train, , drop = FALSE], model = "BayesC"))
    fit <- BGLR::BGLR(y = unname(y[train]), ETA = eta, nIter = nIter,
                      burnIn = burnIn, verbose = verbose, ...)
    u    <- as.numeric(fit$ETA[[1L]]$b)
    gebv <- as.numeric(M %*% u)
  }

  names(gebv) <- d$samples

  # ---- Shape marker effects as a named list of per-trait tables -----------
  # One data.frame per trait (Marker, Allele, Effect); the list name is the
  # trait. This matches the effect_list consumed by run_haplotype_catalogue()
  # and the canonical geno_meta$marker_effects format.
  alt_rows <- data.frame(Marker = amap$Marker, Allele = amap$Alt,
                         Effect = u, stringsAsFactors = FALSE)
  alt_rows <- alt_rows[!is.na(alt_rows$Allele), , drop = FALSE]
  ref_rows <- data.frame(Marker = amap$Marker, Allele = amap$Ref,
                         Effect = 0, stringsAsFactors = FALSE)
  effect_tbl <- rbind(alt_rows, ref_rows)
  effect_tbl <- effect_tbl[order(effect_tbl$Marker, effect_tbl$Allele), ,
                           drop = FALSE]
  rownames(effect_tbl) <- NULL
  marker_effects <- stats::setNames(list(effect_tbl), trait)

  # ---- Shape GEBV phenotype table -----------------------------------------
  gebv_df <- data.frame(Sample = d$samples, stringsAsFactors = FALSE)
  gebv_df[[gebv_colname]] <- gebv[d$samples]

  # ---- Attach to a revalidated GenoMeta copy ------------------------------
  # Merge into any effects already on the GenoMeta so repeated calls accumulate
  # one table per trait; this trait's table replaces any earlier copy of it.
  gm       <- geno_meta
  existing <- .as_effect_list(gm$marker_effects)
  if (is.null(existing)) existing <- list()
  existing[[trait]] <- effect_tbl
  gm$marker_effects <- existing

  # Phenotype: keep the modelled trait values on the GenoMeta beside the GEBV.
  # When the phenotype came from an external source, attach it first so its
  # trait column(s) persist — attach_phenotype() restricts to the shared
  # individuals (genotype samples without a record are NA-filled, external rows
  # for unknown samples dropped). An internal phenotype is already on `gm`. The
  # GEBV column is then merged in, preserving the existing phenotype columns.
  if (!is.null(phenotype))
    gm <- attach_phenotype(gm, phenotype)

  base_ph <- gm$phenotype
  base_ph <- if (is.null(base_ph))
    data.frame(Sample = d$samples, stringsAsFactors = FALSE, check.names = FALSE)
  else
    as.data.frame(base_ph, stringsAsFactors = FALSE, check.names = FALSE)
  base_ph[[gebv_colname]] <- gebv[as.character(base_ph$Sample)]
  gm <- attach_phenotype(gm, base_ph)
  stopifnot(is_geno_meta(gm))

  list(
    marker_effects = marker_effects,
    gebv           = gebv_df,
    model          = fit,
    geno_meta      = gm,
    allele_map     = amap,
    trait          = trait
  )
}


# ---------------------------------------------------------------------------
# Exported: predict whole-genome GEBVs from marker effects
# ---------------------------------------------------------------------------

#' Predict whole-genome GEBVs from marker effects and append them to the phenotype
#'
#' Turns a set of per-marker allele effects into a single whole-genome GEBV per
#' sample and writes it back into the GenoMeta's \code{phenotype} table as a new
#' column, so the observed trait and its predicted GEBV sit side by side and can
#' be compared directly (e.g. with \code{\link{plot_phenotype_bars}}).
#'
#' Unlike \code{\link{fit_genomic_prediction}} — which \emph{estimates} marker
#' effects from a phenotype — this function consumes effects that already exist
#' (supplied externally or attached to the GenoMeta) and only \emph{applies}
#' them. The GEBV is the sum, across the selected block-assigned markers and both
#' haplotypes, of each sample's carried allele effects; it is delegated to
#' \code{\link{compute_profile_localgebv}} (\code{level = "genotype"}) and then
#' summed over blocks, so the same tested effect-join / encoding-validation path
#' is reused.
#'
#' Effect source. An explicit \code{marker_effects} argument always wins;
#' otherwise the effects attached to \code{geno_meta$marker_effects} are used. A
#' \code{message()} states which source was taken. Both are accepted in the
#' canonical \emph{marker-effects} shape — a named list of per-trait
#' \code{(Marker, Allele, Effect)} tables, or a single such table — and are
#' normalised with the same helper the rest of the package uses.
#'
#' Filtering. The \code{markers}, \code{blocks} and \code{chromosomes} arguments
#' restrict which loci contribute to the GEBV (passed straight through to
#' \code{\link{compute_profile_localgebv}}). \strong{Every sample in the GenoMeta
#' is always scored} — there is deliberately no sample filter — so the GEBV
#' column aligns one-to-one with the phenotype rows.
#'
#' @param geno_meta      A valid \code{GenoMeta} list.
#' @param marker_effects Optional marker effects in the canonical shape (named
#'   list of per-trait \code{Marker, Allele, Effect} tables, or a single such
#'   table). \code{NULL} (default) uses \code{geno_meta$marker_effects}.
#' @param trait          Trait whose effect table is applied. \code{NULL}
#'   (default) is allowed only when a single trait is present; with several,
#'   name one. The trait also names the output column.
#' @param gebv_suffix    Suffix for the GEBV column. The column is named
#'   \code{paste0(trait, "_", gebv_suffix)}, so the default \code{"GEBV"} yields
#'   \code{"<trait>_GEBV"}; pass e.g. \code{"GEBV_51b"} for
#'   \code{"<trait>_GEBV_51b"}.
#' @param markers        Optional character vector of marker names to include.
#' @param blocks         Optional character vector of block IDs to include.
#' @param chromosomes    Optional list of chromosome region specs (see
#'   \code{\link{filter_geno_meta}}).
#' @param center        Centring mode applied to the allele effects:
#'   \code{"none"}, \code{"locus_mean"} or \code{"allele_freq"}. See
#'   \code{\link{compute_profile_localgebv}} for what each preserves.
#' @param center_by_freq Logical; mean-centre each allele effect by its
#'   population frequency before summing (forwarded to
#'   \code{\link{compute_profile_localgebv}}). Default \code{FALSE}.
#'
#' @return A copy of \code{geno_meta} with the GEBV column added to its
#'   \code{phenotype} table (existing trait columns preserved) and revalidated
#'   via \code{\link{attach_phenotype}}.
#'
#' @examples
#' \dontrun{
#' meta <- sim_geno_meta(40, 60, block_size = 5L, seed = 1L)
#' ph   <- data.frame(Sample = dimnames(meta$geno_arrays)[[1]],
#'                    Yield  = rnorm(40))
#' meta <- attach_phenotype(meta, ph)
#'
#' # Effects estimated once, then re-applied as a whole-genome GEBV.
#' meta <- fit_genomic_prediction(meta, trait = "Yield")$geno_meta
#' meta <- predict_gebv(meta, trait = "Yield")          # adds "Yield_GEBV"
#'
#' # Compare observed phenotype vs predicted GEBV in a bar chart.
#' plot_phenotype_bars(meta, traits = c("Yield", "Yield_GEBV"))
#'
#' # External effects + a custom suffix.
#' meta <- predict_gebv(meta, marker_effects = my_effect_list,
#'                      trait = "Yield", gebv_suffix = "GEBV_51b")
#' }
#'
#' @seealso \code{\link{fit_genomic_prediction}},
#'   \code{\link{compute_profile_localgebv}}, \code{\link{attach_phenotype}},
#'   \code{\link{plot_phenotype_bars}}
#' @export
predict_gebv <- function(geno_meta,
                         marker_effects = NULL,
                         trait          = NULL,
                         gebv_suffix    = "GEBV",
                         markers        = NULL,
                         blocks         = NULL,
                         chromosomes    = NULL,
                         center         = NULL,
                         center_by_freq = NULL) {
  center <- .resolve_center(center, center_by_freq, "predict_gebv",
                            default = "none")
  if (!is_geno_meta(geno_meta))
    stop("predict_gebv: 'geno_meta' is not a valid GenoMeta.")

  # ---- Resolve the effect source (external overrides the stored effects) ----
  # Both routes go through .as_effect_list() so external input and stored
  # effects are normalised to the SAME canonical named list of per-trait tables.
  if (!is.null(marker_effects)) {
    eff_list <- .as_effect_list(marker_effects)
    message("predict_gebv: using externally supplied marker effects.")
  } else {
    eff_list <- .as_effect_list(geno_meta$marker_effects)
    if (is.null(eff_list))
      stop("predict_gebv: no marker effects supplied and none attached to the ",
           "GenoMeta. Pass 'marker_effects' or run fit_genomic_prediction() ",
           "first.")
    message("predict_gebv: using marker effects stored on the GenoMeta.")
  }
  if (length(eff_list) == 0L)
    stop("predict_gebv: the marker-effect source is empty.")

  # ---- Select one trait's effect table (and name the output column) --------
  if (is.null(trait)) {
    if (length(eff_list) > 1L)
      stop("predict_gebv: marker effects hold multiple traits (",
           paste(names(eff_list), collapse = ", "),
           "); set 'trait' to choose one.")
    trait <- names(eff_list)[1L]
  } else if (!trait %in% names(eff_list)) {
    stop("predict_gebv: trait '", trait, "' not found. Available: ",
         paste(names(eff_list), collapse = ", "), ".")
  }
  effect_df <- eff_list[[trait]]

  # ---- Per-block local GEBV (all samples, requested loci) -> per-sample sum -
  # samples = NULL scores every sample; markers/blocks/chromosomes restrict the
  # loci. Genotype-level columns are samples, so a colSums over blocks gives each
  # sample's whole-genome GEBV.
  block_gebv <- compute_profile_localgebv(
    geno_meta,
    effect_df      = effect_df,
    samples        = NULL,
    markers        = markers,
    blocks         = blocks,
    chromosomes    = chromosomes,
    level          = "genotype",
    center         = center)

  sample_cols <- setdiff(names(block_gebv), "Block")
  if (length(sample_cols) == 0L)
    stop("predict_gebv: no GEBV could be computed — no block-assigned markers ",
         "remain after filtering. Attach block info with attach_block_info() ",
         "or widen the marker/block filter.")
  gebv <- colSums(as.matrix(block_gebv[, sample_cols, with = FALSE]),
                  na.rm = TRUE)

  # ---- Append the GEBV column to the phenotype ([Trait]_[suffix]) ----------
  # Preserve any existing trait columns: build on the stored phenotype (or a
  # bare Sample frame when none is attached) and add the GEBV column aligned by
  # Sample, then revalidate through attach_phenotype().
  gebv_col     <- paste0(trait, "_", gebv_suffix)
  geno_samples <- dimnames(geno_meta$geno_arrays)[[1L]]

  ph <- geno_meta$phenotype
  if (is.null(ph))
    ph <- data.frame(Sample = geno_samples, stringsAsFactors = FALSE,
                     check.names = FALSE)
  else
    ph <- as.data.frame(ph, stringsAsFactors = FALSE, check.names = FALSE)

  ph[[gebv_col]] <- unname(gebv[as.character(ph$Sample)])

  attach_phenotype(geno_meta, ph)
}
