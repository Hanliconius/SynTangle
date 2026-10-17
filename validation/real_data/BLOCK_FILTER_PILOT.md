# Size-filtered search pilot

Uses the completed annotated Lepidoptera fixture and frozen best saved layout. No downloads, BLAST, MCScanX, or GENESPACE reruns.

Five concurrent Pegasus tasks request 180 solver seconds each. Control retains all links for 180 seconds. Treatments remove the smallest 10%, 25%, 50%, or 75% of pair-specific links for 60 seconds, halve that removal fraction for 60 seconds, then restore every link for 60 seconds. Rank by geometric mean of block spans in bp across the two assemblies, with lexical ID tie breaking. This is block-span filtering, not BLAST confidence filtering.

Every chromosome and hard orientation equation remains present. Filtering may expose smaller graph components; component sizes and solver diagnostics are recorded per stage. A reduced graph is rebuilt when evidence changes. No conditional exclusions or filtered optimality claims are transferred to the restored graph. Only complete legal layouts are carried forward.

Every candidate is rescored on the original full fixture, including removed links. The best ordinary-crossing layout cannot regress against the frozen baseline. The sequential trajectory may pass through layouts that improve the filtered objective while worsening full crossings; these do not replace the full-data incumbent. A separate best visual-score layout is retained and can have more ordinary crossings.

The diagnostic span-weighted score sums sqrt(w_i*w_j) over crossing link pairs, where each link weight is the geometric mean of its two block spans divided by their chromosome lengths. This is an exploratory proxy for ribbon prominence, not a certified weighted optimization objective. Absolute crossing counts and full-data bounds remain authoritative for the original objective. Filtered bounds are reported separately and never substituted for a full-data certificate.

Budget parity covers requested solver time, not model construction or diagnostic rescoring; wall time is reported. Jobs use one CPU and 12 GB. No PDF rendering is performed by this first search pilot; saved states support later comparison plots without rerunning preprocessing.

Run from the repository root:

```bash
bash validation/real_data/submit_block_filter.sh local_results/annotated_lep_rjv857
```

Read logs:

```bash
cat logs/st_filter.*.{out,err}
```

Summarize only completed tasks:

```bash
python - <<'PY'
import json
from pathlib import Path
root=Path('local_results/latest_block_filter_run.txt').read_text().strip()
for i in range(5):
    p=Path(root)/f'task_{i}'
    if not (p/'COMPLETE').exists():
        print(i,'INCOMPLETE; inspect task log')
        continue
    d=json.loads((p/'report.json').read_text())
    print('removed=',d['filter_fraction'],'C=',d['final']['crossings'],
          'weighted=',round(d['final']['span_weighted_crossings'],3),
          'L=',d['full_lower_bound'],'seconds=',round(d['wall_seconds'],2),
          'components=',[s['component_sizes'] for s in d['stages']])
PY
```

Validation before publication: Python compilation, shell syntax, and an input-only check of filtering membership, chromosome/constraint preservation, and diagnostic agreement with the canonical crossing scorer. Solver benchmarks must run on Pegasus.
