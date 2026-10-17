# GENESPACE integration from R

The flip helper improves orientations within your existing GENESPACE chromosome order. Full SynTangle instead searches a new chromosome order and orientations from the same evidence. Read the [user guide](walkthrough.md) for the distinction and statistical interpretation.

Install Python SynTangle and the bundled `SynTangleR` package as shown in the main README. Both R functions call the same Python implementation used by the CLI. No homology discovery is repeated and no GENESPACE functions are patched.

## Input

Pass the object returned by `GENESPACE::plot_riparian()`. The adapter reads `plotData$sourceData$chromosomes` and `plotData$sourceData$blocks`; a direct object containing `$sourceData` is also accepted. This adapter targets the GENESPACE 1.3.1 source-data schema. Compatibility checks stop explicitly when required fields are absent.

Preserve all the arguments used to generate the original plot, including `genomeIDs`, `refGenome`, `useRegions`, `useOrder`, and any ordering functions. `useOrder=TRUE` uses gene-rank coordinates; `FALSE` uses base-pair coordinates. The adapter uses GENESPACE's **displayed chromosome extents**, which may differ from full assembly lengths. It does not extend or truncate the source evidence.

```r
library(SynTangleR)
flips <- syntangle_flip_assist(
    rp, useOrder = FALSE, output = "results/flips",
    executable = "/absolute/path/to/SynTangle/.venv/bin/syntangle"
)
print(flips)
flips$stats
flips$pairStats
```

If the original plot already inverted chromosomes:

```r
flips <- syntangle_flip_assist(
    rp, useOrder = FALSE, invertTheseChrs = originalInversions,
    output = "results/flips_from_existing"
)
```

GENESPACE source blocks precede its display inversion transform, so this original inversion list is required to recover the displayed starting state. The adapter cannot reliably infer it from chromosome labels.

## Replot with GENESPACE

For flip assistance, reuse your original `plot_riparian()` arguments and replace `invertTheseChrs` with `flips$invertTheseChrs`. For example, if these were the original plotting arguments:

```r
rp_flipped <- GENESPACE::plot_riparian(
    gsParam = gsParam, useRegions = TRUE, useOrder = FALSE,
    invertTheseChrs = flips$invertTheseChrs
)
```

The returned inversion list is absolute relative to the source data: **replace, do not concatenate with the original list**. SynTangle itself holds chromosome order fixed during flip assistance. GENESPACE runs its ordering code again when replotting; keep all ordering arguments unchanged and confirm the chromosome sequence still matches `flips$layout`, particularly where native ordering ties occur. The comparison PDF from SynTangle draws the exact audited input and output layouts.

## Full optimization

```r
fit <- syntangle_run(
    rp, useOrder = FALSE, output = "results/full", seconds = 600,
    ribbonObjective = "crossing-cap", ribbonSeconds = 300
)
fit$layout
fit$stats
fit$pdf
```

Full optimization can change chromosome order. Passing only `fit$invertTheseChrs` back to GENESPACE **does not apply the optimized order**. Use SynTangle's vector PDF for the optimized drawing; `fit$layout` supplies explicit ranks and orientations for custom renderers. Colors in the SynTangle renderer identify supplied links, not GENESPACE's reference-chromosome palette.

## Returned object

`layout` gives species, chromosome, rank and orientation. `inputLayout` records the original display. `stats` contains crossing counts, link count, lower bound and gap (where available), and status. `pairStats` reports each adjacent species pair. `changedChromosomes` lists relative flips; `invertTheseChrs` is the final absolute inversion list. `audit` and `pdf` give output paths. The complete R object is saved as `syntangle_result.rds`.

## Export only

```r
syntangle_export_genespace(rp, "exported_synteny", useOrder = FALSE)
```

This creates `chromosomes.tsv`, `links.tsv`, and a source description for the CLI. Every region is retained. Missing displayed chromosomes, nonadjacent links, inconsistent coordinates and zero-width intervals cause explicit errors. The public adapter never silently fills in missing GENESPACE orders or removes evidence to make an input pass.

Adapter tests exercise the documented source-data schema and an installed Python engine. An end-to-end run against each new GENESPACE version is needed before claiming compatibility with that version.
