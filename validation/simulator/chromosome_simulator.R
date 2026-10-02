# SynTangle validation simulator
#
# Forward simulator for small, inspectable chromosome-evolution fixtures.
# Biological evolution (fusion/fission/inversion) is deliberately separate from
# display-tangle induction. The solver must never be given the hidden event log.

empty_evolution_log <- function() {
  data.frame(
    event_id = integer(),
    species = character(),
    type = character(),
    input_chromosomes = character(),
    output_chromosomes = character(),
    start = numeric(),
    end = numeric(),
    span_fraction = numeric(),
    detail = character(),
    stringsAsFactors = FALSE
  )
}

simulate_ancestor_genome <- function(
  species = "speciesA",
  n_chrom = 10L,
  genes_per_chrom = 40L,
  chromosome_span = 5e6,
  gene_length_range = c(1000L, 50000L),
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n_chrom >= 1L, genes_per_chrom >= 2L)

  all_genes <- vector("list", n_chrom)
  chroms <- vector("list", n_chrom)
  start_grid <- seq.int(1L, as.integer(chromosome_span), by = 1000L)

  if (length(start_grid) < genes_per_chrom) {
    stop("chromosome_span is too small for genes_per_chrom at 1 kb start-grid spacing")
  }

  for (chr_idx in seq_len(n_chrom)) {
    starts <- sort(sample(start_grid, genes_per_chrom, replace = FALSE))
    lengths <- sample(
      seq.int(as.integer(gene_length_range[1]), as.integer(gene_length_range[2])),
      genes_per_chrom,
      replace = TRUE
    )
    ends <- starts + lengths
    chrom <- paste0("chr", chr_idx)

    all_genes[[chr_idx]] <- data.frame(
      name = paste0("gene", chr_idx, "_", seq_along(starts)),
      start = as.numeric(starts),
      end = as.numeric(ends),
      strand = sample(c("+", "-"), genes_per_chrom, replace = TRUE),
      chrom = chrom,
      species = species,
      stringsAsFactors = FALSE
    )

    chroms[[chr_idx]] <- data.frame(
      chrom = chrom,
      length = max(ends),
      species = species,
      stringsAsFactors = FALSE
    )
  }

  list(
    genes = do.call(rbind, all_genes),
    chroms = do.call(rbind, chroms),
    evolution_log = empty_evolution_log()
  )
}

clone_genome <- function(genome, species, inherit_log = TRUE) {
  out <- genome
  out$genes$species <- species
  out$chroms$species <- species
  if (!inherit_log || is.null(out$evolution_log)) {
    out$evolution_log <- empty_evolution_log()
  }
  out
}

append_evolution_event <- function(
  log,
  event_id,
  species,
  type,
  input_chromosomes,
  output_chromosomes,
  start = NA_real_,
  end = NA_real_,
  span_fraction = NA_real_,
  detail = ""
) {
  row <- data.frame(
    event_id = as.integer(event_id),
    species = as.character(species),
    type = as.character(type),
    input_chromosomes = paste(input_chromosomes, collapse = ";"),
    output_chromosomes = paste(output_chromosomes, collapse = ";"),
    start = as.numeric(start),
    end = as.numeric(end),
    span_fraction = as.numeric(span_fraction),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  rbind(log, row)
}

choose_inversion_window <- function(genes_chr, chr_len, span_range = c(0.4, 0.9)) {
  genes_chr <- genes_chr[order(genes_chr$start, genes_chr$end), , drop = FALSE]
  stopifnot(length(span_range) == 2L, span_range[1] > 0, span_range[2] <= 1,
            span_range[1] < span_range[2])

  # Prefer breakpoints in true intergenic gaps so no simulated gene is cut.
  internal_cuts <- numeric()
  if (nrow(genes_chr) >= 2L) {
    gap_ok <- genes_chr$end[-nrow(genes_chr)] < genes_chr$start[-1] - 1
    left_end <- genes_chr$end[-nrow(genes_chr)][gap_ok]
    right_start <- genes_chr$start[-1][gap_ok]
    if (length(left_end)) {
      internal_cuts <- floor((left_end + right_start) / 2)
    }
  }

  cuts <- sort(unique(c(1, internal_cuts, chr_len)))
  if (length(cuts) < 2L) stop("No legal inversion breakpoints available")

  candidates <- list()
  k <- 1L
  for (i in seq_len(length(cuts) - 1L)) {
    for (j in (i + 1L):length(cuts)) {
      inv_start <- cuts[i]
      inv_end <- cuts[j]
      frac <- (inv_end - inv_start + 1) / chr_len
      if (frac < span_range[1] || frac > span_range[2]) next
      contained <- genes_chr$start >= inv_start & genes_chr$end <= inv_end
      if (!any(contained)) next
      candidates[[k]] <- c(start = inv_start, end = inv_end, span_fraction = frac)
      k <- k + 1L
    }
  }

  if (!length(candidates)) {
    stop("No intergenic inversion window satisfies the requested span range")
  }

  candidates[[sample.int(length(candidates), 1L)]]
}

mutate_genome_structure <- function(
  genome,
  n_mutations = 3L,
  event_plan = NULL,
  allowed_events = c("fusion", "fission", "inversion"),
  inversion_span = c(0.4, 0.9),
  fusion_gap = 10000,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  genes <- genome$genes
  chroms <- genome$chroms
  log <- if (is.null(genome$evolution_log)) empty_evolution_log() else genome$evolution_log
  species <- unique(genes$species)
  if (length(species) != 1L) stop("A genome passed to mutate_genome_structure must contain exactly one species")

  available_chroms <- unique(chroms$chrom)
  used_chroms <- character()

  if (!is.null(event_plan)) {
    if (!all(event_plan %in% c("fusion", "fission", "inversion"))) {
      stop("event_plan contains an unknown event type")
    }
    n_mutations <- length(event_plan)
  } else {
    n_mutations <- as.integer(n_mutations)
  }

  if (n_mutations < 0L) stop("n_mutations must be >= 0")

  next_event_id <- if (nrow(log)) max(log$event_id) + 1L else 1L

  for (i in seq_len(n_mutations)) {
    possible <- setdiff(available_chroms, used_chroms)
    if (!length(possible)) break

    if (!is.null(event_plan)) {
      action <- event_plan[i]
    } else {
      feasible <- intersect(allowed_events, c(
        if (length(possible) >= 2L) "fusion",
        if (length(possible) >= 1L) "fission",
        if (length(possible) >= 1L) "inversion"
      ))
      if (!length(feasible)) break
      action <- sample(feasible, 1L)
    }

    event_id <- next_event_id
    next_event_id <- next_event_id + 1L

    if (action == "fusion") {
      if (length(possible) < 2L) stop("Requested fusion but fewer than two unused chromosomes remain")
      chr_pair <- sample(possible, 2L, replace = FALSE)
      used_chroms <- c(used_chroms, chr_pair)

      g1 <- genes[genes$chrom == chr_pair[1], , drop = FALSE]
      g2 <- genes[genes$chrom == chr_pair[2], , drop = FALSE]
      g1 <- g1[order(g1$start), , drop = FALSE]
      g2 <- g2[order(g2$start), , drop = FALSE]

      shift <- max(g1$end) + fusion_gap - min(g2$start) + 1
      g2$start <- g2$start + shift
      g2$end <- g2$end + shift

      new_chrom <- paste0("fusion_", chr_pair[1], "_", chr_pair[2], "_", event_id)
      g1$chrom <- new_chrom
      g2$chrom <- new_chrom
      fused <- rbind(g1, g2)

      genes <- genes[!(genes$chrom %in% chr_pair), , drop = FALSE]
      genes <- rbind(genes, fused)
      chroms <- chroms[!(chroms$chrom %in% chr_pair), , drop = FALSE]
      chroms <- rbind(
        chroms,
        data.frame(chrom = new_chrom, length = max(fused$end), species = species,
                   stringsAsFactors = FALSE)
      )

      log <- append_evolution_event(
        log, event_id, species, "fusion", chr_pair, new_chrom,
        detail = paste0("gap=", fusion_gap)
      )
    }

    if (action == "fission") {
      splittable <- possible[vapply(possible, function(chr) sum(genes$chrom == chr) >= 2L, logical(1))]
      if (!length(splittable)) stop("Requested fission but no unused chromosome has at least two genes")
      chr_split <- sample(splittable, 1L)
      used_chroms <- c(used_chroms, chr_split)

      g <- genes[genes$chrom == chr_split, , drop = FALSE]
      g <- g[order(g$start), , drop = FALSE]
      split_index <- floor(nrow(g) / 2)
      f1 <- g[seq_len(split_index), , drop = FALSE]
      f2 <- g[(split_index + 1L):nrow(g), , drop = FALSE]

      f1_name <- paste0(chr_split, "_f1_", event_id)
      f2_name <- paste0(chr_split, "_f2_", event_id)
      f1$chrom <- f1_name
      f2$chrom <- f2_name

      # Historical workflow required the second fission product to begin at 1.
      f2_offset <- min(f2$start) - 1
      f2$start <- f2$start - f2_offset
      f2$end <- f2$end - f2_offset

      genes <- genes[genes$chrom != chr_split, , drop = FALSE]
      genes <- rbind(genes, f1, f2)
      chroms <- chroms[chroms$chrom != chr_split, , drop = FALSE]
      chroms <- rbind(
        chroms,
        data.frame(chrom = f1_name, length = max(f1$end), species = species,
                   stringsAsFactors = FALSE),
        data.frame(chrom = f2_name, length = max(f2$end), species = species,
                   stringsAsFactors = FALSE)
      )

      log <- append_evolution_event(
        log, event_id, species, "fission", chr_split, c(f1_name, f2_name),
        start = max(f1$end),
        detail = paste0("f2_coordinate_offset=", f2_offset)
      )
    }

    if (action == "inversion") {
      invertible <- possible[vapply(possible, function(chr) sum(genes$chrom == chr) >= 2L, logical(1))]
      if (!length(invertible)) stop("Requested inversion but no unused chromosome has at least two genes")
      chr_inv <- sample(invertible, 1L)
      used_chroms <- c(used_chroms, chr_inv)

      g <- genes[genes$chrom == chr_inv, , drop = FALSE]
      chr_len <- chroms$length[match(chr_inv, chroms$chrom)]
      window <- choose_inversion_window(g, chr_len, inversion_span)
      inv_start <- unname(window["start"])
      inv_end <- unname(window["end"])
      span_fraction <- unname(window["span_fraction"])

      idx <- genes$chrom == chr_inv & genes$start >= inv_start & genes$end <= inv_end
      if (!any(idx)) stop("Selected inversion window contains no genes")

      old_start <- genes$start[idx]
      old_end <- genes$end[idx]
      genes$start[idx] <- inv_start + (inv_end - old_end)
      genes$end[idx] <- inv_start + (inv_end - old_start)
      genes$strand[idx] <- ifelse(genes$strand[idx] == "+", "-", "+")

      log <- append_evolution_event(
        log, event_id, species, "inversion", chr_inv, chr_inv,
        start = inv_start,
        end = inv_end,
        span_fraction = span_fraction,
        detail = paste0("n_genes=", sum(idx))
      )
    }
  }

  genes <- genes[order(genes$chrom, genes$start, genes$end), , drop = FALSE]
  chroms <- chroms[order(chroms$chrom), c("chrom", "length", "species"), drop = FALSE]
  rownames(genes) <- NULL
  rownames(chroms) <- NULL
  rownames(log) <- NULL

  list(genes = genes, chroms = chroms, evolution_log = log)
}

generate_synteny_blocks <- function(genome1, genome2) {
  g1 <- genome1$genes
  g2 <- genome2$genes
  shared <- intersect(g1$name, g2$name)
  if (!length(shared)) {
    return(data.frame(
      homology_id = character(), chr1 = character(), start1 = numeric(), end1 = numeric(),
      chr2 = character(), start2 = numeric(), end2 = numeric(), strand = character(),
      species1 = character(), species2 = character(), stringsAsFactors = FALSE
    ))
  }

  a <- g1[match(shared, g1$name), c("name", "chrom", "start", "end", "strand", "species")]
  b <- g2[match(shared, g2$name), c("name", "chrom", "start", "end", "strand", "species")]

  data.frame(
    homology_id = shared,
    chr1 = a$chrom,
    start1 = a$start,
    end1 = a$end,
    chr2 = b$chrom,
    start2 = b$start,
    end2 = b$end,
    strand = ifelse(a$strand == b$strand, "+", "-"),
    species1 = a$species,
    species2 = b$species,
    stringsAsFactors = FALSE
  )
}

write_syntenyplotter_blocks <- function(blocks, path) {
  cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "strand", "species1", "species2")
  write.table(blocks[, cols, drop = FALSE], path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

write_genome_tables <- function(genomes, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(genomes)) {
    genome <- genomes[[nm]]
    write.table(genome$genes, file.path(directory, paste0(nm, "_genes.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(genome$chroms, file.path(directory, paste0(nm, "_chromosomes.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(genome$evolution_log, file.path(directory, paste0(nm, "_evolution_log.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
