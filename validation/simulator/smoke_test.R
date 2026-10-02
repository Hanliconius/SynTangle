args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "validation/simulator/smoke_test.R"
root <- normalizePath(file.path(dirname(this_file), "../.."), mustWork = TRUE)
source(file.path(root, "validation/simulator/chromosome_simulator.R"))
source(file.path(root, "validation/simulator/tangle_induction.R"))
source(file.path(root, "validation/simulator/visualize_simulation.R"))

A <- simulate_ancestor_genome("speciesA", n_chrom = 8, genes_per_chrom = 30, seed = 101)
B0 <- clone_genome(A, "speciesB")
B <- mutate_genome_structure(B0, event_plan = c("fission"), seed = 102)
C0 <- clone_genome(B, "speciesC")
C <- mutate_genome_structure(C0, event_plan = c("fusion", "inversion"), seed = 103)
genomes <- list(speciesA = A, speciesB = B, speciesC = C)

# Homology identity is preserved through all structural mutations.
stopifnot(setequal(A$genes$name, B$genes$name))
stopifnot(setequal(B$genes$name, C$genes$name))
stopifnot(!anyDuplicated(A$genes$name), !anyDuplicated(B$genes$name), !anyDuplicated(C$genes$name))

# The fission coordinate reset is explicit and testable.
f2 <- B$chroms$chrom[grepl("_f2_", B$chroms$chrom)]
stopifnot(length(f2) == 1L)
stopifnot(min(B$genes$start[B$genes$chrom == f2]) == 1)

# Planned events were actually recorded.
stopifnot(identical(tail(B$evolution_log$type, 1), "fission"))
stopifnot(identical(tail(C$evolution_log$type, 2), c("fusion", "inversion")))
stopifnot(any(C$evolution_log$span_fraction >= 0.4 & C$evolution_log$span_fraction <= 0.9, na.rm = TRUE))

# Inversion changes orientation for at least one inherited anchor.
merged <- merge(B$genes[, c("name", "strand")], C$genes[, c("name", "strand")],
                by = "name", suffixes = c("_B", "_C"))
stopifnot(any(merged$strand_B != merged$strand_C))

AB <- generate_synteny_blocks(A, B)
BC <- generate_synteny_blocks(B, C)
stopifnot(nrow(AB) == nrow(A$genes))
stopifnot(nrow(BC) == nrow(B$genes))

native <- make_native_display_state(genomes)
tangled <- induce_display_tangle(genomes, mode = "strong", seed = 104)
stopifnot(nrow(tangled$display_state) ==
            sum(vapply(genomes, function(g) nrow(g$chroms), integer(1))))
stopifnot(nrow(tangled$tangle_log) == nrow(tangled$display_state))
for (sp in unique(native$species)) {
  stopifnot(setequal(native$chrom[native$species == sp],
                     tangled$display_state$chrom[tangled$display_state$species == sp]))
}

pdf_path <- file.path(tempdir(), "syntangle_simulator_audit.pdf")
write_simulation_audit_pdf(genomes, native, tangled$display_state, pdf_path)
stopifnot(file.exists(pdf_path), file.info(pdf_path)$size > 0)

cat("SynTangle validation simulator smoke test: PASS\n")
