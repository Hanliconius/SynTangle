# Completed mirror/policy comparison

The supplied collector stdout contains 210 unique completed tasks: 126 runs at
150 seconds and 84 at 600 seconds. 168 report matching global bounds. All 84
largest-case 600-second runs reach C=lower=8377; none of the 42 largest-case
150-second runs closes its gap. All 84 small/medium runs report matching bounds
at C=341 / 1258.

The separate unanchored strict-improvement audit reports infeasibility verified
for every one of the nine presentations, including the three largest. This
closes the earlier 120-second audit's unresolved largest-case results under the
600-second-per-positive-component audit allowance. It is a separate model check
using the same numerical backend, not an independent rational certificate.
Raw component records and layouts have not been inspected from this upload.

## What improved

This experiment improves efficiency; it does not expand the biological sample
or add proof coverage beyond the earlier 600-second optimization runs. The
successful strict audits add corroboration for those reported optima. Nine
presentations represent three biology rungs with mild/strong/random displays.

Public-start largest-case pooled median wall times at 600 seconds are:

| Method | Median wall seconds |
|---|---:|
| Direct control | 274.71 |
| Three-species control | 280.91 |
| Adaptive control | 299.46 |
| Mirror direct | 271.81 |
| Mirror three-species | 227.35 |
| Mirror adaptive | 239.56 |
| Mirror selective | 220.04 |

These are six runs per method (three presentations, two repeats), not a general
ranking. Selective's pooled median is about 20% lower than direct control.
Matched paired **global-stage** times give more modest evidence for the mirror
fixing itself: mirror direct is faster in 7/12 pairs, mirror three in 10/12,
mirror adaptive in 10/12 across both starts. Median paired control/mirror
global-time ratios are 1.041, 1.084 and 1.062 respectively. Some pairs regress;
ratios range roughly 0.65–3.12. Detailed wall times are absent from the printed
table, and stage times overlap, so no matched wall-time ratio is manufactured.

At 150 seconds, selective preserves nearly the same largest-case public-start
scores as adaptive control with shorter neighborhood work:

| Presentation | Direct C | Adaptive C | Selective C | Adaptive neighborhood s | Selective neighborhood s |
|---|---:|---:|---:|---:|---:|
| Mild | 74306 | 16103 | 16103 | 31.80 | 19.69 |
| Strong | 74306 | 22696 | 22741 | 33.70 | 21.41 |
| Random | 74306 | 13550 | 13550 | 37.69 | 19.65 |

Selective leaves about 69–82% fewer crossings than direct control here,
but that large benefit was already present with adaptive
neighborhoods. It is not a new step change caused by mirror symmetry or the
policy. Smaller cases are already fast with direct search; neighborhood preludes
often add overhead. Selective skips that work on smaller models.

Mirror variants can use more memory: peak reported RSS is 789.3 MiB, versus about
443.8 MiB across controls in this report. Both are well below the same 16G
request. Increased requested memory is not the source of improvement.

## Evidence and remaining decisions

- [Exact uploaded collector text](collector_stdout.txt).
- [210 parsed records](records.json) and [28 overview rows](overview.json).
- [Reported audit outcomes](audit_status.json) and [provenance/limitations](metadata.json).

The upload does not print tested commit, old job IDs or result-directory identity.
The concurrently submitted jobs 138406327–138406330 must not be assigned to this
earlier report. Collector globbing can combine multiple runs; this upload contains
one complete 210-task report.

Retain exact mirror reduction and selective scheduling as supported experimental
candidates. Do not promote heuristic thresholds to biological scaling laws.
A broad default still needs independent biology seeds, real data, matched
GENESPACE freedoms/budgets, legal visual reconstruction checks and raw proof
record review. No production default changes accompany this archive.
