# Reproducing the code checks

Use Python 3.12 in a fresh environment, then from the repository root:

```bash
python -m pip install -r validation/requirements-publication.txt
bash validation/run_publication_checks.sh /absolute/path/to/check-results
```

The runner records the exact commit and installed versions, runs both suites, and exits nonzero if either fails. It does not download scientific datasets or run the full analysis pipeline. The solver requirements and matplotlib==3.11.2 are pinned to the reported successful environment. This is a direct-dependency specification, not a complete transitive lockfile. Preserve installed_versions.txt with each validation record.

On 14 October 2026 the Pegasus publication audit's validation checkout passed the unit suite after solver dependency installation, and passed all 43 real-data fixture tests after matplotlib installation (job 138484501). Follow-up job 138496210 confirmed that checkout as commit 7e5eb5e161df6d000769d563eb4b0e72c97fb73b and reported numpy 2.3.5, scipy 1.17.0, cvxpy 1.9.3, scs 3.3.1, highspy 1.15.1, and matplotlib 3.11.2. Those results apply to that checkout, not automatically to newer commits or the recovery files.

## Code recovery disposition

All 9,648 files under source_0/source_1 in recovery branch pegasus-publication-recovery-20261014T183358Z are within local_results. This includes installed third-party tool trees and run-specific files. Do not merge that recovery tree into maintained source. Preserve reproducibility through maintained source, dependency declarations, accession/reference manifests, and commands. Individual generated scripts should only be promoted after comparison with their generator and evidence that independent maintained logic is missing.

Passing fixture tests establishes the tested behavior, not full published-data reproduction. Whole-chromosome movement/flip rules remain unchanged.
