write_validation_bundle <- function(
  genomes,
  display_state,
  directory,
  fixture_id = "simulated_validation_case",
  title = "Forward-simulated SynTangle validation case",
  purpose = "Exercise the solver on extant chromosomes generated from hidden structural history"
) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  species_ids <- names(genomes)
  if (is.null(species_ids) || any(!nzchar(species_ids))) {
    stop("genomes must be a named list in intended species-layer order")
  }

  species_table <- data.frame(
    species_id = species_ids,
    species_rank = seq_along(species_ids),
    stringsAsFactors = FALSE
  )

  chromosome_rows <- list()
  occurrence_rows <- list()
  cidx <- 1L
  oidx <- 1L

  for (species_name in species_ids) {
    genome <- genomes[[species_name]]
    if (length(unique(genome$genes$species)) != 1L) {
      stop("Each genome must contain exactly one species")
    }

    state <- display_state[display_state$species == species_name, , drop = FALSE]
    if (!nrow(state)) {
      stop(paste("No display state rows for", species_name))
    }

    for (i in seq_len(nrow(genome$chroms))) {
      chrom <- genome$chroms$chrom[i]
      rank <- state$display_rank[match(chrom, state$chrom)]
      if (is.na(rank)) {
        stop(paste("Missing display rank for", species_name, chrom))
      }
      chromosome_rows[[cidx]] <- data.frame(
        species_id = species_name,
        chromosome_id = chrom,
        length = genome$chroms$length[i],
        display_rank = as.integer(rank),
        stringsAsFactors = FALSE
      )
      cidx <- cidx + 1L
    }

    for (i in seq_len(nrow(genome$genes))) {
      gene <- genome$genes[i, , drop = FALSE]
      occurrence_rows[[oidx]] <- data.frame(
        occurrence_id = paste(species_name, gene$name, sep = ":"),
        homology_id = gene$name,
        species_id = species_name,
        chromosome_id = gene$chrom,
        start = gene$start,
        end = gene$end,
        strand = gene$strand,
        stringsAsFactors = FALSE
      )
      oidx <- oidx + 1L
    }
  }

  chromosomes <- do.call(rbind, chromosome_rows)
  occurrences <- do.call(rbind, occurrence_rows)
  chromosomes <- chromosomes[
    order(match(chromosomes$species_id, species_ids), chromosomes$display_rank),
    ,
    drop = FALSE
  ]
  occurrences <- occurrences[
    order(
      match(occurrences$species_id, species_ids),
      occurrences$chromosome_id,
      occurrences$start,
      occurrences$end
    ),
    ,
    drop = FALSE
  ]

  metadata <- data.frame(
    key = c("fixture_id", "title", "purpose"),
    value = c(fixture_id, title, purpose),
    stringsAsFactors = FALSE
  )

  write.table(
    species_table,
    file.path(directory, "species.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    chromosomes,
    file.path(directory, "chromosomes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    occurrences,
    file.path(directory, "occurrences.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    metadata,
    file.path(directory, "metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # Hidden-truth outputs are colocated for benchmarking convenience but are
  # deliberately not consumed by the Python bundle importer.
  evolution_logs <- lapply(genomes, function(genome) genome$evolution_log)
  evolution_logs <- evolution_logs[vapply(evolution_logs, nrow, integer(1)) > 0L]
  if (length(evolution_logs)) {
    hidden_evolution <- unique(do.call(rbind, evolution_logs))
  } else {
    hidden_evolution <- empty_evolution_log()
  }

  write.table(
    hidden_evolution,
    file.path(directory, "hidden_evolution_log.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  input_state <- display_state[
    order(match(display_state$species, species_ids), display_state$display_rank),
    ,
    drop = FALSE
  ]
  write.table(
    input_state,
    file.path(directory, "input_display_state.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  invisible(directory)
}

write_hidden_tangle_log <- function(tangle_log, directory) {
  write.table(
    tangle_log,
    file.path(directory, "hidden_tangle_log.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(directory)
}
