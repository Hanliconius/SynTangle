# Release preparation

This is a development readiness record, not a release announcement.

## Implemented

- Explicit JSON evidence / complete TSV layout contract for standalone refinement.
- Flips-only and full adaptive modes with legal-state validation and rescoring.
- Four-panel biological comparisons retain evidence when native orders omit chromosomes.
- Supplemented baselines are labelled and record every appended chromosome.
- Tests cover invalid imports, coupled orientation constraints, failed baselines,
  partial-order completion and directed endpoint preservation.
- Installed-wheel CI exercises the public entry point outside the source tree.

## Required before tagging

- Verify updated GENESPACE R wrapper on Pegasus and archive its installed function.
- Independently replay certificates used for mathematical optimality claims.
- Review matched benchmark inputs, endpoint geometry and each source reconstruction.
- Complete version-specific R-object import/export adapter, or explicitly omit it
  from the first release and support only the documented TSV interface.
- Confirm third-party/source-data redistribution terms; the package MIT license is included.
- Remove development-only artifacts from built distributions; keep historical evidence
  accessible without requiring users to download genomes or archived result ZIPs.
- Run a clean install and all supported CLI modes in CI before tagging a version.

## Current biological evidence

Suzukii and Rhodotorula completed runs are recorded in Pegasus output and saved
comparison archives. Rhodotorula uses minimap2 block reconstruction and supplemented
GENESPACE ordering, not an end-to-end GENESPACE gene-level run. Ficus is a partial
set of pairwise reconstructions. Faba is blocked on a defensible fixed representation
of split chromosome 1. Do not present these scopes as exact original figure replication.

The previous README is preserved under docs/development rather than silently
rewriting historical algorithm and benchmark statements.

## Local checks during preparation

- 147 Python core tests passed with solver and experimental test dependencies.
- Targeted adapter/endpoint tests passed, including normal partial completion and
  genuine external failure handling.
- Wheel built successfully; a clean virtual environment installed it without optional
  dependencies and ran CLI help, fixture analysis and flips-only refinement outside
  the checkout.
- Installed R wrapper still requires the Pegasus adapter-only check. No release tag
  or claim of a complete GENESPACE R-object integration has been made.
