# helper-ld-toy.R — planted-LD toy panels for the LD-resolver tests
# (test-selectiontools-blocks.R, test-hapselect-blocks.R). testthat sources
# every helper-*.R before the tests.

# Toy panel with a KNOWN block structure: groups of identical markers (r^2 = 1
# within a group) separated by independent draws (r^2 ~ 1/n between groups), on
# chromosomes named 1H..3H. Lines are inbred (A/A or T/T).
st_toy <- function(sizes = c(4L, 1L, 3L, 1L, 5L, 2L, 1L, 4L, 3L, 1L, 2L),
                   chrom = rep(c("1H", "2H", "3H"), c(5L, 3L, 3L)),
                   n = 40L, seed = 11L, noise_in = NA_integer_) {
  set.seed(seed)
  cells <- NULL; mk <- NULL; ch <- NULL; pos <- NULL
  p <- 0
  for (g in seq_along(sizes)) {
    base <- sample(c("A", "T"), n, replace = TRUE)
    for (k in seq_len(sizes[g])) {
      v <- base
      if (!is.na(noise_in) && g == noise_in && k == ceiling(sizes[g] / 2))
        v <- sample(c("A", "T"), n, replace = TRUE)          # breaks the LD run
      cells <- cbind(cells, paste0(v, "/", v))
      mk  <- c(mk, sprintf("g%02d_m%d", g, k))
      ch  <- c(ch, chrom[g])
      p   <- p + 1000
      pos <- c(pos, p)
    }
    p <- p + 1e6
  }
  dimnames(cells) <- list(sprintf("L%02d", seq_len(n)), mk)
  map <- data.frame(Chrom = ch, Marker = mk, Pos = pos, stringsAsFactors = FALSE)
  suppressMessages(build_geno_meta(cells, map, split_char = "/",
                                   encoding = "nucleotide"))
}
planted_group <- function(marker) sub("_m[0-9]+$", "", marker)
partition <- function(bm) {
  s <- split(as.character(bm$Marker), as.character(bm$Block))
  sort(vapply(s, function(x) paste(sort(x), collapse = "+"), ""))
}
