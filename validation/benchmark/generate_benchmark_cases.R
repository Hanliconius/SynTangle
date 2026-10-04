args <- commandArgs(trailingOnly = TRUE)
valid_profiles <- c(
  "smoke", "full", "paired-smoke", "paired", "scale", "stress", "factorial"
)
if (length(args) < 1L || length(args) > 2L) {
  stop(
    paste0(
      "Usage: Rscript validation/benchmark/generate_benchmark_cases.R ",
      "OUTPUT_ROOT [", paste(valid_profiles, collapse = "|"), "]"
    )
  )
}
output_root <- args[[1]]
profile <- if (length(args) == 2L) args[[2]] else "smoke"
if (!profile %in% valid_profiles) {
  stop(paste("profile must be one of", paste(valid_profiles, collapse = ", ")))
}

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
this_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  "validation/benchmark/generate_benchmark_cases.R"
}
repo_root <- normalizePath(
  file.path(dirname(this_file), "../.."),
  mustWork = TRUE
)

source(file.path(repo_root, "validation/simulator/chromosome_simulator.R"))
source(file.path(repo_root, "validation/simulator/tangle_induction.R"))
source(file.path(repo_root, "validation/simulator/export_validation_bundle.R"))

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

event_plan_for_branch <- function(branch_index, intensity) {
  if (intensity == "none") {
    plans <- list(character())
  } else if (intensity == "low") {
    plans <- list(
      c("fission"),
      c("fusion"),
      c("inversion")
    )
  } else if (intensity == "medium") {
    plans <- list(
      c("fission", "inversion"),
      c("fusion", "inversion"),
      c("fission", "fusion")
    )
  } else if (intensity == "high") {
    # Five independent whole-chromosome structural events per lineage step.
    # Two fusions + two fissions keep expected chromosome count roughly stable.
    plans <- list(
      c("fusion", "fission", "inversion", "fusion", "fission"),
      c("fission", "fusion", "inversion", "fission", "fusion"),
      c("inversion", "fusion", "fission", "fusion", "fission")
    )
  } else if (intensity == "very_high") {
    # Eight events per lineage step: three fusion/fission pairs plus two
    # inversions. At least eleven distinct chromosomes are required.
    plans <- list(
      c(
        "fusion", "fission", "inversion", "fusion",
        "fission", "inversion", "fusion", "fission"
      ),
      c(
        "fission", "fusion", "inversion", "fission",
        "fusion", "inversion", "fission", "fusion"
      ),
      c(
        "inversion", "fusion", "fission", "fusion",
        "inversion", "fission", "fusion", "fission"
      )
    )
  } else if (intensity == "extreme") {
    # Twelve events per lineage step: four fusion/fission pairs plus four
    # inversions. This deliberately creates much more cumulative structural
    # history while keeping chromosome number approximately stable.
    plans <- list(
      c(
        "fusion", "fission", "inversion",
        "fusion", "fission", "inversion",
        "fusion", "fission", "inversion",
        "fusion", "fission", "inversion"
      ),
      c(
        "fission", "fusion", "inversion",
        "fission", "fusion", "inversion",
        "fission", "fusion", "inversion",
        "fission", "fusion", "inversion"
      ),
      c(
        "inversion", "fusion", "fission",
        "inversion", "fusion", "fission",
        "inversion", "fusion", "fission",
        "inversion", "fusion", "fission"
      )
    )
  } else {
    stop(paste("Unknown event intensity:", intensity))
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
  seed,
  biology_id = case_id,
  biology_seed = seed,
  tangle_seed = seed + 9000L,
  lineage_model = "sequential",
  benchmark_axis = "legacy",
  benchmark_level = ""
) {
  species_ids <- paste0("species", LETTERS[seq_len(n_species)])

  if (!lineage_model %in% c("sequential", "independent")) {
    stop("lineage_model must be sequential or independent")
  }

  ancestor_species <- if (lineage_model == "independent") {
    "__ancestor__"
  } else {
    species_ids[[1]]
  }

  ancestor <- simulate_ancestor_genome(
    species = ancestor_species,
    n_chrom = n_chrom,
    genes_per_chrom = genes_per_chrom,
    chromosome_span = 1500000,
    seed = biology_seed
  )

  genomes <- list()

  if (lineage_model == "sequential") {
    genomes[[species_ids[[1]]]] <- ancestor

    for (branch_index in 2:n_species) {
      previous <- genomes[[species_ids[[branch_index - 1L]]]]
      descendant <- clone_genome(previous, species_ids[[branch_index]])
      descendant <- mutate_genome_structure(
        descendant,
        event_plan = event_plan_for_branch(branch_index, intensity),
        seed = biology_seed + branch_index * 100L
      )
      genomes[[species_ids[[branch_index]]]] <- descendant
    }
  } else {
    # Orthogonal scaling benchmarks use independent extant descendants from a
    # common hidden ancestor. This prevents increasing species count from also
    # increasing cumulative rearrangement depth.
    for (branch_index in seq_len(n_species)) {
      descendant <- clone_genome(
        ancestor,
        species_ids[[branch_index]],
        inherit_log = FALSE
      )
      descendant <- mutate_genome_structure(
        descendant,
        event_plan = event_plan_for_branch(branch_index + 1L, intensity),
        seed = biology_seed + branch_index * 100L
      )
      genomes[[species_ids[[branch_index]]]] <- descendant
    }
  }

  native <- make_native_display_state(genomes)
  tangled <- induce_display_tangle(
    genomes,
    mode = tangle_mode,
    seed = tangle_seed,
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
    purpose = paste0(
      "Measure legal untangling performance on hidden-truth forward ",
      "simulations"
    )
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
    biology_id = biology_id,
    species_count = n_species,
    ancestor_chromosomes = n_chrom,
    genes_per_chromosome = genes_per_chrom,
    tangle_mode = tangle_mode,
    event_intensity = intensity,
    events_per_branch = length(event_plan_for_branch(2L, intensity)),
    seed = seed,
    biology_seed = biology_seed,
    tangle_seed = tangle_seed,
    lineage_model = lineage_model,
    benchmark_axis = benchmark_axis,
    benchmark_level = benchmark_level,
    stringsAsFactors = FALSE
  )
}

spec <- function(
  case_id,
  n_species,
  n_chrom,
  genes_per_chrom,
  tangle_mode,
  intensity,
  seed,
  biology_id = case_id,
  biology_seed = seed,
  tangle_seed = seed + 9000L,
  lineage_model = "sequential",
  benchmark_axis = "legacy",
  benchmark_level = ""
) {
  list(
    case_id,
    n_species,
    n_chrom,
    genes_per_chrom,
    tangle_mode,
    intensity,
    seed,
    biology_id,
    biology_seed,
    tangle_seed,
    lineage_model,
    benchmark_axis,
    benchmark_level
  )
}

paired_specs <- function(
  species_values,
  chromosome_values,
  prefix,
  genes_per_chrom = 12L,
  seed_base = 5000L
) {
  specs <- list()
  k <- 1L
  biology_index <- 1L
  mode_offsets <- c(mild = 9101L, strong = 9201L, random = 9301L)

  for (n_species in species_values) {
    for (n_chrom in chromosome_values) {
      intensity <- if (n_chrom >= 8L) "medium" else "low"
      biology_id <- sprintf(
        "%s_%dsp_%dchr",
        prefix,
        n_species,
        n_chrom
      )
      biology_seed <- seed_base + biology_index * 101L

      for (tangle_mode in c("mild", "strong", "random")) {
        case_id <- sprintf(
          "%s_%s",
          biology_id,
          tangle_mode
        )
        # Keep the optimizer seed fixed within a biological triplet so
        # differences among modes reflect presentation state rather than a
        # different stochastic restart sequence.
        solver_seed <- biology_seed
        tangle_seed <- biology_seed + mode_offsets[[tangle_mode]]

        specs[[k]] <- spec(
          case_id,
          n_species,
          n_chrom,
          genes_per_chrom,
          tangle_mode,
          intensity,
          solver_seed,
          biology_id = biology_id,
          biology_seed = biology_seed,
          tangle_seed = tangle_seed
        )
        k <- k + 1L
      }

      biology_index <- biology_index + 1L
    }
  }

  specs
}

if (profile == "smoke") {
  specs <- list(
    spec(
      "smoke_3sp_4chr_mild",
      3L, 4L, 8L, "mild", "low", 1301L
    ),
    spec(
      "smoke_3sp_5chr_strong",
      3L, 5L, 8L, "strong", "medium", 1302L
    ),
    spec(
      "smoke_4sp_5chr_random",
      4L, 5L, 6L, "random", "medium", 1303L
    )
  )
} else if (profile == "full") {
  # Historical Stage 13/15 benchmark. Seeds remain case-specific so this
  # profile stays reproducible and directly comparable with existing results.
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
        case_seed <- 2000L + k * 37L
        specs[[k]] <- spec(
          case_id,
          n_species,
          n_chrom,
          12L,
          tangle_mode,
          intensity,
          case_seed
        )
        k <- k + 1L
      }
    }
  }
} else if (profile == "paired-smoke") {
  specs <- paired_specs(
    species_values = c(3L),
    chromosome_values = c(5L),
    prefix = "paired_smoke",
    genes_per_chrom = 8L,
    seed_base = 4000L
  )
} else if (profile == "paired") {
  # Controlled presentation-tangle experiment: each mild/strong/random triplet
  # is generated from exactly the same evolved extant genomes.
  specs <- paired_specs(
    species_values = c(3L, 4L, 5L),
    chromosome_values = c(6L, 8L, 10L),
    prefix = "paired",
    genes_per_chrom = 12L,
    seed_base = 5000L
  )
} else if (profile == "scale") {
  # First scaling probe beyond the original 10-chromosome grid. Species count
  # is held at four while chromosome count increases; each scale point is a
  # controlled mild/strong/random triplet.
  specs <- paired_specs(
    species_values = c(4L),
    chromosome_values = c(12L, 15L, 20L, 30L),
    prefix = "scale",
    genes_per_chrom = 12L,
    seed_base = 7000L
  )
} else if (profile == "factorial") {
  # Biologically grounded orthogonal benchmark. Extant species are independent
  # descendants of one hidden ancestor, so changing species number does not
  # also add sequential rearrangement depth. Presentation mode is fixed at
  # random: paired presentation invariance is tested elsewhere.
  #
  # The design deliberately centers the chromosome axis on 16 (a useful
  # general eukaryotic reference from the Genome Observatory discussion) and
  # 31 (a common Lepidoptera-like karyotype), while retaining 40 only as a
  # high-end boundary probe.
  factorial_cases <- list(
    # Species-count axis: chromosome count and rearrangement burden fixed.
    list("factor_species_4sp_16chr_1ev", 4L, 16L, "low", 9101L, "species", "4"),
    list("factor_species_8sp_16chr_1ev", 8L, 16L, "low", 9102L, "species", "8"),
    list("factor_species_12sp_16chr_1ev", 12L, 16L, "low", 9103L, "species", "12"),
    list("factor_species_20sp_16chr_1ev", 20L, 16L, "low", 9104L, "species", "20"),

    # Chromosome-count axis: species number and rearrangement burden fixed.
    # The 8sp/16chr/1ev point above is the shared baseline and is not repeated.
    list("factor_chrom_8sp_8chr_1ev", 8L, 8L, "low", 9201L, "chromosomes", "8"),
    list("factor_chrom_8sp_31chr_1ev", 8L, 31L, "low", 9202L, "chromosomes", "31"),
    list("factor_chrom_8sp_40chr_1ev", 8L, 40L, "low", 9203L, "chromosomes", "40"),

    # Rearrangement axis: species and chromosome count fixed.
    # Again, the 1-event point is the shared 8sp/16chr baseline.
    list("factor_events_8sp_16chr_0ev", 8L, 16L, "none", 9301L, "rearrangements", "0"),
    list("factor_events_8sp_16chr_2ev", 8L, 16L, "medium", 9302L, "rearrangements", "2"),
    list("factor_events_8sp_16chr_5ev", 8L, 16L, "high", 9303L, "rearrangements", "5"),

    # Lepidoptera-like conserved-karyotype checks.
    list("factor_lep_8sp_31chr_0ev", 8L, 31L, "none", 9401L, "lepidoptera_like", "8sp_0ev"),
    list("factor_lep_12sp_31chr_0ev", 12L, 31L, "none", 9402L, "lepidoptera_like", "12sp_0ev"),
    list("factor_lep_12sp_31chr_1ev", 12L, 31L, "low", 9403L, "lepidoptera_like", "12sp_1ev")
  )

  specs <- lapply(
    factorial_cases,
    function(item) {
      case_id <- item[[1]]
      n_species <- item[[2]]
      n_chrom <- item[[3]]
      intensity <- item[[4]]
      biology_seed <- item[[5]]
      benchmark_axis <- item[[6]]
      benchmark_level <- item[[7]]
      spec(
        case_id,
        n_species,
        n_chrom,
        12L,
        "random",
        intensity,
        biology_seed,
        biology_id = case_id,
        biology_seed = biology_seed,
        tangle_seed = biology_seed + 9301L,
        lineage_model = "independent",
        benchmark_axis = benchmark_axis,
        benchmark_level = benchmark_level
      )
    }
  )
} else if (profile == "stress") {
  # Coupled stress ladder: species count, chromosome count, and cumulative
  # structural-event density all increase together. Each rung is represented
  # by a mild/strong/random presentation triplet of the same evolved biology.
  stress_rungs <- list(
    list(6L, 20L, 16L, "high", 5L, 8101L),
    list(8L, 30L, 18L, "very_high", 8L, 8201L),
    list(10L, 40L, 20L, "extreme", 12L, 8301L)
  )

  specs <- list()
  k <- 1L
  mode_offsets <- c(mild = 9101L, strong = 9201L, random = 9301L)

  for (rung in stress_rungs) {
    n_species <- rung[[1]]
    n_chrom <- rung[[2]]
    genes_per_chrom <- rung[[3]]
    intensity <- rung[[4]]
    expected_events <- rung[[5]]
    biology_seed <- rung[[6]]
    biology_id <- sprintf(
      "stress_%dsp_%dchr_%dev",
      n_species,
      n_chrom,
      expected_events
    )

    if (length(event_plan_for_branch(2L, intensity)) != expected_events) {
      stop("Stress rung event count does not match intensity definition")
    }

    for (tangle_mode in c("mild", "strong", "random")) {
      case_id <- sprintf("%s_%s", biology_id, tangle_mode)
      specs[[k]] <- spec(
        case_id,
        n_species,
        n_chrom,
        genes_per_chrom,
        tangle_mode,
        intensity,
        biology_seed,
        biology_id = biology_id,
        biology_seed = biology_seed,
        tangle_seed = biology_seed + mode_offsets[[tangle_mode]]
      )
      k <- k + 1L
    }
  }
}

manifest_rows <- lapply(
  specs,
  function(item) {
    do.call(
      make_case,
      setNames(
        item,
        c(
          "case_id",
          "n_species",
          "n_chrom",
          "genes_per_chrom",
          "tangle_mode",
          "intensity",
          "seed",
          "biology_id",
          "biology_seed",
          "tangle_seed",
          "lineage_model",
          "benchmark_axis",
          "benchmark_level"
        )
      )
    )
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

cat(
  "Generated",
  nrow(manifest),
  profile,
  "benchmark cases in",
  output_root,
  "\n"
)
