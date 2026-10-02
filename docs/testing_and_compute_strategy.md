# Validation and compute strategy

SynTangle will use three validation tiers.

## Tier A — synthetic executable fixtures

**Location:** Git repository / GitHub Actions  
**Purpose:** correctness and invariants  
**Scale:** tiny

Synthetic cases have known answers by construction. They test representation, decomposition, orientation propagation, admissible order constraints, crossing calculations, and the prohibition on subchromosomal reordering.

These tests should run on every pull request and never require HPC resources.

## Tier B — reduced real-data smoke tests

**Location:** Git repository if genuinely small and redistributable; otherwise referenced externally  
**Purpose:** confirm that the canonical data model and parsers behave sensibly on actual synteny output  
**Scale:** small

Once the core representation/decomposition code exists, add one or more curated real examples. These should be reduced enough that unexpected behavior can still be inspected manually.

A real dataset is not a substitute for Tier A because the biological ground truth is generally unknown.

## Tier C — full real-data / scaling benchmarks

**Location:** Pegasus working directory  
**Purpose:** performance, scaling, realistic optimization, large graph behavior  
**Scale:** potentially large

Full synteny matrices, whole-genome block tables, large parameter sweeps, spectral calculations on large unresolved kernels, optimizer restarts, and scaling benchmarks should run on Pegasus rather than in GitHub Actions.

Git should contain:

- source code;
- configuration;
- Slurm launch scripts;
- compact manifests/checksums;
- summarized benchmark results;
- small diagnostic outputs when useful.

Git should not contain large working datasets, large intermediate graphs, optimizer caches, or bulk benchmark output.

## Development boundary

The default progression is:

```text
prove correctness on synthetic fixtures
        ↓
smoke-test on a small real dataset
        ↓
profile locally / in CI
        ↓
move computationally substantial work to Pegasus
        ↓
commit only reproducible scripts + compact results
```

No algorithm should be judged correct solely because it produces a visually plausible answer on real data.
