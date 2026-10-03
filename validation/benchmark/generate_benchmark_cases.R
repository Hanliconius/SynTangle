args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: Rscript validation/benchmark/generate_benchmark_cases.R OUTPUT_ROOT [smoke|full]")
}
output_root <- args[[1]]
profile <- if (length(args) == 2L) args[[2]] else "smoke"
if (!profile %in% c("smoke", "full")) {
  stop("profile must be smoke or full")
}

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
this_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  "validation/benchmark/generate_benchmark_cases.R"
}
repo_root <- normalizePath(file.path(dirname(this_file), "../.."), mustWork = TRUE)

source(file.path(repo_root, "validation/simulator/chromosome_simulator.R"))
source(file.path(repo_root, "validation/simulator/tangle_induction.R"))
source(file.path(repo_root, "validation/simulator/export_validation_bundle.R"))

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

event_plan_for_branch <- function(branch_index, intensity) {
  if (intensity == "low") {
    plans <- list(
      c("fission"),
      c("fusion"),
      c("inversion")
    )
  } else {
    plans <- list(
      c("fission", "inversion"),
      c("fusion", "inversion"),
      c("fission", "fusion")
    )
  }
  plans[[((branch_index - 2L) %% length(plans)) + 1L]]
}

make_case <- function(
  case_id,
  n_species,
  n_chrom,
  genes_per_chrom,
  tangle_mode,
  intensity,
  seed
) {
  species_ids <- paste0("species", LETTERS[seq_len(n_species)])

  ancestor <- simulate_ancestor_genome(
    species = species_ids[[1]],
    n_chrom = n_chrom,
    genes_per_chrom = genes_per_chrom,
    chromosome_span = 1500000,
    seed = seed
  )

  genomes <- list()
  genomes[[species_ids[[1]]]] <- ancestor

  for (branch_index in 2:n_species) {
    previous <- genomes[[species_ids[[branch_index - 1L]]]]
    descendant <- clone_genome(previous, species_ids[[branch_index]])
    descendant <- mutate_genome_structure(
      descendant,
      event_plan = event_plan_for_branch(branch_index, intensity),
      seed = seed + branch_index * 100L
    )
    genomes[[species_ids[[branch_index]]]] <- descendant
  }

  native <- make_native_display_state(genomes)
  tangled <- induce_display_tangle(
    genomes,
    mode = tangle_mode,
    seed = seed + 9000L,
    flip_probability = switch(
      tangle_mode,
      mild = 0.15,
      strong = 0.35,
      random = 0.50
    )
  )

  case_dir <- file.path(output_root, case_id)
  write_validation_bundle(
    genomes = genomes,
    display_state = tangled$display_state,
    directory = case_dir,
    fixture_id = case_id,
    title = paste("SynTangle simulation benchmark", case_id),
    purpose = "Measure legal untangling performance on hidden-truth forward simulations"
  )
  write_hidden_tangle_log(tangled$tangle_log, case_dir)

  write.table(
    native,
    file.path(case_dir, "hidden_native_display_state.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  data.frame(
    case_id = case_id,
    case_dir = case_id,
    species_count = n_species,
    ancestor_chromosomes = n_chrom,
    genes_per_chromosome = genes_per_chrom,
    tangle_mode = tangle_mode,
    event_intensity = intensity,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

if (profile == "smoke") {
  specs <- list(
    list("smoke_3sp_4chr_mild", 3L, 4L, 8L, "mild", "low", 1301L),
    list("smoke_3sp_5chr_strong", 3L, 5L, 8L, "strong", "medium", 1302L),
    list("smoke_4sp_5chr_random", 4L, 5L, 6L, "random", "medium", 1303L)
  )
} else {
  specs <- list()
  k <- 1L
  for (n_species in c(3L, 4L, 5L)) {
    for (n_chrom in c(6L, 8L, 10L)) {
      for (tangle_mode in c("mild", "strong", "random")) {
        intensity <- if (n_chrom >= 8L) "medium" else "low"
        case_id <- sprintf(
          "full_%dsp_%dchr_%s",
          n_species,
          n_chrom,
          tangle_mode
        )
        specs[[k]] <- list(
          case_id,
          n_species,
          n_chrom,
          12L,
          tangle_mode,
          intensity,
          2000L + k * 37L
        )
        k <- k + 1L
      }
    }
  }
}

manifest_rows <- lapply(
  specs,
  function(spec) {
    do.call(make_case, setNames(
      spec,
      c(
        "case_id",
        "n_species",
        "n_chrom",
        "genes_per_chrom",
        "tangle_mode",
        "intensity",
        "seed"
      )
    ))
  }
)
manifest <- do.call(rbind, manifest_rows)

write.table(
  manifest,
  file.path(output_root, "benchmark_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Generated", nrow(manifest), profile, "benchmark cases in", output_root, "\n")
