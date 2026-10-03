# Forward chromosome simulator

This directory contains small forward simulations used to validate SynTangle against known biological truth.

The simulator is **validation machinery, not part of the SynTangle solver**. The solver should see only extant chromosomes/homology and an initial display state. It must not see the hidden evolutionary event log or the deliberate tangle log.

## Architecture

```text
ancestor genome
    ↓
fusion / fission / inversion
    ↓
extant descendant genomes
    ↓
hidden evolution log

extant descendant genomes
    ↓
whole-chromosome display permutation / reversal only
    ↓
tangled presentation supplied to SynTangle
    ↓
hidden tangle log
```

The two histories are intentionally separate. Biological simulation may change chromosome membership and within-chromosome coordinates according to explicit structural mutations. Tangle induction may change only whole-chromosome display rank and display orientation. Both rank and orientation are part of the public initial layout supplied to SynTangle; the operations that produced them remain hidden in the tangle log.

## Files

- `chromosome_simulator.R` — cleaned forward simulator, event log, synteny-anchor export.
- `tangle_induction.R` — seeded presentation-only permutation/reversal.
- `visualize_simulation.R` — lightweight base-R PDF audit comparing native and deliberately tangled presentations.
- `smoke_test.R` — tiny invariant test for CI.
- `export_validation_bundle.R` — writes public extant-state TSVs plus separate hidden-truth logs.
- `generate_integration_case.R` — deterministic fission/fusion/inversion case used by CI.
- `check_integration_case.py` — imports only the public bundle, runs SynTangle, then reveals hidden history for benchmark checks.
- `archive/reconstructed_2025_simulator.R` — provenance copy of the recovered 2025 logic; this is historical reference, not production validation code.

## Historical provenance

The simulator was reconstructed from the June 2025 development conversation. The recovered implementation used stable gene identities across descendants and applied chromosome fusion, fission and inversion while retaining a mutation log. The old code also became the source of the chromosome-similarity / bipartite-ordering experiments that ultimately motivated SynTangle.

The cleaned implementation deliberately fixes several historical rough edges:

1. random mutation counts are scalar rather than relying on `sample(2:4)`;
2. the second fission product is recentered using one saved offset for both `start` and `end`;
3. inversion breakpoints are preferentially placed in intergenic gaps, so a simulated gene is not silently cut at a breakpoint;
4. fusion spacing is an explicit gap between the two participating chromosomes;
5. event plans can be supplied to make validation fixtures deterministic.

None of these changes alter SynTangle's normative rules. They only make the forward validation data easier to reason about and reproduce.

## Minimal example

```r
source("validation/simulator/chromosome_simulator.R")
source("validation/simulator/tangle_induction.R")
source("validation/simulator/visualize_simulation.R")

A <- simulate_ancestor_genome("speciesA", seed = 1)
B <- mutate_genome_structure(clone_genome(A, "speciesB"),
                             event_plan = c("fission"), seed = 2)
C <- mutate_genome_structure(clone_genome(B, "speciesC"),
                             event_plan = c("fusion", "inversion"), seed = 3)

genomes <- list(speciesA = A, speciesB = B, speciesC = C)
native <- make_native_display_state(genomes)
tangled <- induce_display_tangle(genomes, mode = "strong", seed = 4)

write_simulation_audit_pdf(
  genomes,
  native,
  tangled$display_state,
  "simulation_audit.pdf"
)
```

The standard development principle is to inspect the audit plot as well as machine assertions. A graph-theoretically valid synthetic case is not useful if its simulated chromosome biology is nonsensical.
