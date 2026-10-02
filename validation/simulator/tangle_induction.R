# SynTangle display-tangle induction
#
# These functions deliberately change only the PRESENTATION state. They never
# edit genomic coordinates, homology, or chromosome membership.

make_native_display_state <- function(genomes) {
  rows <- list()
  k <- 1L
  for (nm in names(genomes)) {
    chroms <- genomes[[nm]]$chroms
    for (i in seq_len(nrow(chroms))) {
      rows[[k]] <- data.frame(
        species = unique(chroms$species)[1],
        chrom = chroms$chrom[i],
        source_rank = i,
        display_rank = i,
        orientation = 1L,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

swap_adjacent <- function(x, n_swaps) {
  if (length(x) < 2L || n_swaps <= 0L) return(x)
  for (i in seq_len(n_swaps)) {
    j <- sample.int(length(x) - 1L, 1L)
    x[c(j, j + 1L)] <- x[c(j + 1L, j)]
  }
  x
}

induce_display_tangle <- function(
  genomes,
  mode = c("mild", "strong", "random"),
  seed = NULL,
  mild_swaps = 2L,
  flip_probability = NULL
) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(flip_probability)) {
    flip_probability <- switch(mode, mild = 0.15, strong = 0.35, random = 0.50)
  }
  if (flip_probability < 0 || flip_probability > 1) stop("flip_probability must be between 0 and 1")

  native <- make_native_display_state(genomes)
  out <- native

  for (sp in unique(native$species)) {
    idx <- which(native$species == sp)
    chroms <- native$chrom[idx]
    order_chr <- switch(
      mode,
      mild = swap_adjacent(chroms, mild_swaps),
      strong = sample(chroms, length(chroms), replace = FALSE),
      random = sample(chroms, length(chroms), replace = FALSE)
    )

    out$display_rank[idx] <- match(native$chrom[idx], order_chr)
    out$orientation[idx] <- ifelse(runif(length(idx)) < flip_probability, -1L, 1L)
  }

  out <- out[order(out$species, out$display_rank), , drop = FALSE]
  rownames(out) <- NULL

  log <- merge(
    native[, c("species", "chrom", "display_rank", "orientation")],
    out[, c("species", "chrom", "display_rank", "orientation")],
    by = c("species", "chrom"),
    suffixes = c("_before", "_after"),
    sort = FALSE
  )
  log$flipped <- log$orientation_before != log$orientation_after
  log$mode <- mode
  log <- log[, c("species", "chrom", "display_rank_before", "display_rank_after", "flipped", "mode")]

  list(display_state = out, tangle_log = log)
}
