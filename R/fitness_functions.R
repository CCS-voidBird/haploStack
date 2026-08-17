# =============================================================================
# fitness_functions.R
#
# create_fitness_index()  — aggregate marker effects into a composite fitness
#   score across multiple traits with configurable weights.
# update_fitness_index()  — update an existing fitness index with new data.
# write_fitness_file()    — serialise a fitness index to disk.
# run_fitness_pipeline()  — high-level wrapper for the full fitness workflow.
# =============================================================================

#' Create a composite fitness index from multi-trait marker effects
#'
#' Combines per-marker allele effects across multiple traits into a single
#' composite "fitness" effect table using trait weights. The result has the
#' same structure as a single-trait effect file (Marker, Allele, Effect) and
#' is used as input to the stacking simulation.
#'
#' @param effect_list    Named list of effect data.frames, one per trait. Each
#'   must have at least three columns: Marker, Allele, Effect (read by
#'   position).
#' @param trait_weights  Named numeric vector of weights, one per trait. Names
#'   must match names in \code{effect_list}.
#'
#'   \strong{Trait direction is the user's responsibility.} The composite simply
#'   sums \code{weight * effect} across traits, so a \emph{lower-is-better} trait
#'   (e.g. a disease score) must be given a \strong{negative} weight to flip its
#'   direction — otherwise selection drives it the wrong way. For a
#'   value-magnitude summary, use the absolute weights to judge each trait's
#'   share. Example: \code{c(Yield = 1, DiseaseScore = -0.5)} rewards yield and
#'   penalises disease.
#'
#'   Weights are used \strong{as given} by default, so the composite is the plain
#'   \code{sum(weight * effect)} in the traits' own units. Pass
#'   \code{normalise_weights = TRUE} to rescale them first.
#' @param allele_dict    Character vector of valid allele states (default:
#'   \code{DEFAULT_ALLELE_DICT}). Used to build the complete marker × allele
#'   registry.
#' @param normalise_weights Logical. \code{FALSE} (default) uses the weights
#'   exactly as supplied, so the composite effects — and therefore any GEBV
#'   computed from them — stay in raw units. \code{TRUE} divides the weights by
#'   \code{sum(abs(weights))} so the absolute weights sum to 1.0 (preserving sign
#'   for direction), which rescales the composite by a positive constant: useful
#'   for comparing indices built from different weight vectors, and harmless to
#'   selection, since a positive constant cannot change any ranking.
#' @param fill_missing   Numeric value used for (Marker, Allele) combinations
#'   absent from a trait's effect table. Default \code{0}.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{\code{Marker}}{Marker ID (character).}
#'   \item{\code{Allele}}{Allele state (character).}
#'   \item{\code{Effect}}{Composite weighted effect (numeric).}
#'   \item{\code{<trait>}}{One column per input trait showing its contribution.}
#' }
#'
#' @examples
#' eff1 <- data.frame(Marker = c("M1","M1","M2","M2"),
#'                    Allele = c("A","T","A","T"),
#'                    Effect = c(0.5, -0.5, 0.2, -0.2))
#' eff2 <- data.frame(Marker = c("M1","M1","M2","M2"),
#'                    Allele = c("A","T","A","T"),
#'                    Effect = c(0.1, -0.1, 0.8, -0.8))
#' fitness <- create_fitness_index(
#'   effect_list   = list(Yield = eff1, Protein = eff2),
#'   trait_weights = c(Yield = 0.7, Protein = 0.3)
#' )
#' print(fitness)
#'
#' @seealso \code{\link{update_fitness_index}}, \code{\link{write_fitness_file}},
#'   \code{\link{run_fitness_pipeline}}
#'
#' @export
create_fitness_index <- function(effect_list,
                                  trait_weights,
                                  allele_dict       = DEFAULT_ALLELE_DICT,
                                  normalise_weights = FALSE,
                                  fill_missing      = 0) {

  # Validate inputs
  if (length(effect_list) == 0L)
    stop("create_fitness_index: 'effect_list' must contain at least one trait.")
  trait_names <- names(effect_list)
  if (is.null(trait_names))
    stop("create_fitness_index: 'effect_list' must be named.")
  if (!all(trait_names %in% names(trait_weights)))
    stop("create_fitness_index: 'trait_weights' must have names matching all traits.")

  # Weights are used AS GIVEN by default (0.5.4): the composite is the plain
  # sum(weight * effect), so its effects — and any GEBV derived from them — stay
  # in raw units. With normalise_weights = TRUE they are divided by the sum of
  # their ABSOLUTE values, so that a signed weight vector (negative = "lower is
  # better", see @param trait_weights) still sums to unit magnitude and each
  # trait's contribution stays balanced. Dividing by sum(w) would blow up / flip
  # when weights sum to ~0 (e.g. +1 and -1).
  w <- trait_weights[trait_names]
  if (any(is.na(w))) stop("create_fitness_index: some trait weights are NA.")
  if (normalise_weights) {
    denom <- sum(abs(w))
    if (denom == 0) stop("create_fitness_index: 'trait_weights' are all zero.")
    w <- w / denom
  }

  # Standardise each effect table (read by position)
  eff_dfs <- lapply(seq_along(effect_list), function(i) {
    df <- as.data.frame(effect_list[[i]], stringsAsFactors = FALSE)
    if (ncol(df) < 3L)
      stop(sprintf("Effect table for '%s' needs >= 3 columns.", trait_names[i]))
    colnames(df)[1:3] <- c("Marker", "Allele", "Effect")
    df$Marker <- as.character(df$Marker)
    df$Allele <- as.character(df$Allele)
    df$Effect <- as.numeric(df$Effect)
    df[, c("Marker","Allele","Effect"), drop = FALSE]
  })
  names(eff_dfs) <- trait_names

  # Build complete marker × allele registry:
  # all unique markers that appear in any effect table × all allele_dict values
  all_markers <- unique(unlist(lapply(eff_dfs, function(df) df$Marker)))
  all_pairs <- CJ(Marker = all_markers, Allele = allele_dict)

  # Compute weighted composite
  result <- all_pairs

  trait_cols <- list()
  for (trait in trait_names) {
    eff_dt <- as.data.table(eff_dfs[[trait]])
    joined <- merge(result, eff_dt, by = c("Marker","Allele"), all.x = TRUE)
    joined$Effect[is.na(joined$Effect)] <- fill_missing
    trait_cols[[trait]] <- joined$Effect * w[[trait]]
  }

  result$Effect <- Reduce("+", trait_cols)
  for (trait in trait_names) {
    result[[trait]] <- trait_cols[[trait]] / w[[trait]]  # un-weighted contribution
  }

  # Sort by Marker, Allele
  result <- result[order(result$Marker, result$Allele), ]
  rownames(result) <- NULL

  as.data.table(result)
}


#' Update a fitness index with new effect data
#'
#' Re-computes the composite fitness effect for a given set of traits using
#' updated effect estimates, while preserving the marker-allele registry from
#' the original fitness index.
#'
#' @param fitness_index  data.table returned by \code{\link{create_fitness_index}}.
#' @param updated_effects Named list of updated effect data.frames (same
#'   structure as \code{effect_list} in \code{\link{create_fitness_index}}).
#' @param trait_weights   Named numeric vector of weights (same trait names).
#' @param normalise_weights Logical. Default \code{FALSE} (weights used as
#'   given; see \code{\link{create_fitness_index}}).
#'
#' @return Updated \code{data.table} with the same structure as
#'   \code{fitness_index}.
#'
#' @export
update_fitness_index <- function(fitness_index,
                                  updated_effects,
                                  trait_weights,
                                  normalise_weights = FALSE) {

  if (!is.data.table(fitness_index))
    stop("update_fitness_index: 'fitness_index' must be a data.table.")
  if (!all(c("Marker","Allele","Effect") %in% names(fitness_index)))
    stop("update_fitness_index: 'fitness_index' must have Marker, Allele, Effect columns.")

  # Rebuild using the restricted marker × allele set from the existing index
  base_registry <- fitness_index[, .(Marker, Allele)]
  allele_dict   <- unique(fitness_index$Allele)

  create_fitness_index(
    effect_list       = updated_effects,
    trait_weights     = trait_weights,
    allele_dict       = allele_dict,
    normalise_weights = normalise_weights
  )
}


#' Write a fitness index to disk
#'
#' Serialises a fitness index data.table to a space-delimited text file
#' suitable for use as an effect file in the stacking simulation.
#'
#' @param fitness_index data.table from \code{\link{create_fitness_index}}.
#' @param path          Output file path.
#' @param cols          Character vector of columns to write. Default writes
#'   only \code{c("Marker","Allele","Effect")}.
#' @param sep           Column separator (default \code{" "}).
#'
#' @return Invisible \code{NULL}.
#'
#' @export
write_fitness_file <- function(fitness_index,
                                path,
                                cols = c("Marker","Allele","Effect"),
                                sep  = " ") {
  cols_exist <- intersect(cols, names(fitness_index))
  if (length(cols_exist) == 0L)
    stop("write_fitness_file: none of the requested columns found in fitness_index.")

  fwrite(fitness_index[, cols_exist, with = FALSE],
         file = path, sep = sep, quote = FALSE, row.names = FALSE)
  message("Fitness index written to: ", path)
  invisible(NULL)
}


#' Run the full fitness index pipeline
#'
#' High-level wrapper that creates, optionally visualises, and saves a fitness
#' index from a list of trait effect files.
#'
#' @param effect_list    Named list of effect data.frames (Marker, Allele, Effect).
#' @param trait_weights  Named numeric vector of trait weights.
#' @param output_file    Path to write the fitness index. When \code{NULL}
#'   (default) the index is returned without writing.
#' @param allele_dict    Valid allele states (default \code{DEFAULT_ALLELE_DICT}).
#' @param normalise_weights Logical. Default \code{FALSE} (weights used as
#'   given; see \code{\link{create_fitness_index}}).
#' @param plot_landscape Logical. When \code{TRUE}, print a barplot of
#'   marker-level fitness effects. Default \code{FALSE}.
#'
#' @return The fitness index \code{data.table}, invisibly.
#'
#' @export
run_fitness_pipeline <- function(effect_list,
                                  trait_weights,
                                  output_file       = NULL,
                                  allele_dict       = DEFAULT_ALLELE_DICT,
                                  normalise_weights = FALSE,
                                  plot_landscape    = FALSE) {

  message("Building fitness index for ", length(effect_list), " trait(s)...")
  fi <- create_fitness_index(
    effect_list       = effect_list,
    trait_weights     = trait_weights,
    allele_dict       = allele_dict,
    normalise_weights = normalise_weights
  )
  message("  Fitness index: ", nrow(fi), " marker-allele rows.")

  if (plot_landscape) {
    p <- ggplot(fi, aes(x = reorder(paste0(Marker,".",Allele), Effect),
                        y = Effect, fill = Effect > 0)) +
      geom_col() +
      scale_fill_manual(values = c("TRUE" = "#2ca25f", "FALSE" = "#de2d26"),
                        guide  = "none") +
      labs(title = "Fitness Landscape — Marker Effects",
           x = "Marker.Allele", y = "Composite Effect") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
    print(p)
  }

  if (!is.null(output_file)) {
    write_fitness_file(fi, output_file)
  }

  invisible(fi)
}
