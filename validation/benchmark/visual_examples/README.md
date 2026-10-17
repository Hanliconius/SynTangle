# Saved-layout visual summaries

Submit on Pegasus:

```bash
bash validation/benchmark/submit_visual_summary.sh
```

The nano job snapshots code and the saved-layout ZIP. It produces a five-page
vector PDF, PNG previews and `layout_verification.json`. No optimization,
GENESPACE run or proof audit is repeated. No image-generation service is used.
It uses the existing `syntangle_test` environment as the base for an isolated
plotting environment with matplotlib 3.10.8.

Every saved state is legally validated and canonically re-scored. The figures
compare three random stress presentations: input, native GENESPACE pilot output,
our flip assistance, and SynTangle short/long results where available. A separate
synthetic zero example demonstrates 1188 avoidable crossings. Small/medium cases
already close their gaps under the short allowance; no missing long run is invented.

Within-chromosome coordinates, links, colours and scales are retained within each
case. Colours repeat beyond 20 homology groups. SynTangle starts from the frozen
public state independently of GENESPACE; these panels are not a sequential
pipeline. All rows may move in SynTangle; GENESPACE uses reference-based variants.
Difficulty labels concern these fixtures, not a general graph classification.
Historical runtimes and pooled medians have different timing scopes; inclusive
stage timings overlap and are not summed.

Override the archive with `SCT_VIS_ZIP=/absolute/path/to/bundle.zip`, preserving
the bundled directory structure. To render an extracted archive directly:

```bash
PYTHONPATH=src python validation/benchmark/visual_summary.py \
  --data-root /path/to/extracted_bundle --output /path/to/figures
```
