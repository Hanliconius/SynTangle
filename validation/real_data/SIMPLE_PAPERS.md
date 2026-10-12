# Four smaller published synteny cases

Selected panels: Ficus Fig. 4; faba/pea/Medicago Fig. 2b (all five rows); Rhodotorula linear ribbon panel (scope to verify); suzukii/melanogaster Fig. 4. Circular comparisons are excluded.

## Stage 1: source acquisition

`bash validation/real_data/submit_simple_papers.sh intake` submits four independent, unthrottled CPU tasks and a collector. It saves pinned deposit assets, checks published MD5 hashes, retrieves supplementary figures/tables and repository inventories, and records unresolved source requirements. Existing assets are reused only after checksum validation. Partial downloads never enter the input cache. No raw reads are required; faba's 3.4 GB genome is deliberately excluded because deposited proteins and GFF suffice for a protein-based reconstruction.

The initial run is **source acquisition**, not a completed plot comparison. Comparator assembly versions, author links, and orientations remain to be established. Rhodotorula metadata are collected for all deposited candidates; no arbitrary strain subset is selected. The faba code archive contains general analysis scripts rather than a published synteny plotting script. A reconstruction must be labelled explicitly.

## Stage 2: reviewed conversion

Each case's `prepared/` directory must contain `fixture.json`, `provenance.json`, `input_layout.json`, and `source_review.json`. Fixture schema is the standard SynTangle schema. Full chromosome lengths, signed occurrences, source link identities and displayed exclusions must be preserved. Input layout uses the standard saved-layout encoding.

`source_review.json` must record status PASS, assembly_versions_verified=true, published_layout_verified=true, source_kind (published_evidence or reconstructed_homology), and sha256 for the three input files. A source review should also record original file hashes, selection rules and supporting figure evidence. No hashes alone establish scientific correctness. These gates prevent unreviewed layouts being presented as the published panel.

Set provenance input_panel_title/input_layout_note to distinguish exact evidence replication from a reconstructed homology display in recovered paper order. The shared renderer now supports explicit validated `input_layout.json`; older cases retain their original behavior.

## Stage 3: comparison

`bash validation/real_data/submit_simple_papers.sh compare RUN` checks all four reviewed inputs before submitting twelve concurrent, single-core comparisons at 30, 120 and 600 seconds of SynTangle solver time. Full workflow wall time and peak RSS are recorded separately; these are not assumed equal to solver budgets. GENESPACE ordering uses the existing block-proxy workflow and our flip assistance, all rescored on identical evidence. Endpoint preservation checks are mandatory. Collection reports missing runs rather than claiming success. Independent proof certification remains on the separate exact-certification-experiment branch and is not implied by these comparison scripts.

Output: per-case source_audit.json, comparisons/BUDGET/comparison.pdf, saved layouts and audits, per-task resource measurements, summary.json and SynTangle_simple_papers.zip. Data preparation or matching is performed on Pegasus, not during local code tests.
