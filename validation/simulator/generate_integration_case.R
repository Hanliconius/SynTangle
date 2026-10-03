args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript validation/simulator/generate_integration_case.R OUTPUT_DIR")
}
output_dir <- args[[1]]

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
this_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  "validation/simulator/generate_integration_case.R"
}
root <- normalizePath(file.path(dirname(this_file), "../.."), mustWork = TRUE)

source(file.path(root, "validation/simulator/chromosome_simulator.R"))
source(file.path(root, "validation/simulator/tangle_induction.R"))
source(file.path(root, "validation/simulator/export_validation_bundle.R"))

A <- simulate_ancestor_genome(
  "speciesA",
  n_chrom = 4,
  genes_per_chrom = 8,
  chromosome_span = 800000,
  seed = 801
)

B <- mutate_genome_structure(
  clone_genome(A, "speciesB"),
  event_plan = c("fission"),
  seed = 802
)

C <- mutate_genome_structure(
  clone_genome(B, "speciesC"),
  event_plan = c("fusion", "inversion"),
  seed = 803
)

genomes <- list(speciesA = A, speciesB = B, speciesC = C)
tangled <- induce_display_tangle(
  genomes,
  mode = "strong",
  seed = 804,
  flip_probability = 0.35
)

write_validation_bundle(
  genomes = genomes,
  display_state = tangled$display_state,
  directory = output_dir,
  fixture_id = "simulator_ci_fission_fusion_inversion"
)
write_hidden_tangle_log(tangled$tangle_log, output_dir)

cat(output_dir, "\n")
