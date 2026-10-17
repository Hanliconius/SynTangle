# Experimental JCVI and minimap2 importers

**Status: implemented and tested on small synthetic records; not yet validated on real JCVI or minimap2 datasets.** Do not treat converted figures as verified reproductions of an upstream plot. Inspect the import audit and endpoint correspondence before using results in a manuscript.

Both importers produce `chromosomes.tsv`, `links.tsv`, `excluded.tsv` and `import_audit.json`. They require chromosome-size files for both species; these are two-column `chromosome length` tables (FASTA `.fai` files are also accepted). Row order sets the initial chromosome order. All chromosomes in the size files are retained, including those without links. The two species remain in query-then-target order.

## JCVI / MCscan

The supported input is the standard **six-column `.simple` block file**, plus the corresponding BED gene coordinates. Raw `.anchors`, `--rich` and `--coords` variants are not accepted directly. Generate a simple file with JCVI's `compara.synteny simple` command using the BED files from that analysis.

```bash
syntangle import-jcvi \
  --simple query.target.simple \
  --query-bed query.bed --target-bed target.bed \
  --query-sizes query.fa.fai --target-sizes target.fa.fai \
  --query-species query --target-species target \
  --output imported_jcvi
```

Each simple row has four endpoint gene IDs, a numeric block score, and `+` or `-` relative orientation. The importer forms the genomic envelope of the endpoint BED genes on each side. It applies the explicit relative orientation to the target interval. Gene transcription strands do not override the block orientation. BED coordinates are interpreted as zero-based, half-open intervals.

Every block row becomes a separate link, even if two rows describe the same coordinates. Missing genes, duplicate BED gene IDs, interchromosomal block endpoints or inconsistent endpoint ordering stop conversion. JCVI highlight/color prefixes are accepted but their styling is not carried into SynTangle. This importer reproduces block intervals, not JCVI's exact graphical geometry or colors.

Format reference: [JCVI SimpleFile and simple exporter](https://github.com/tanghaibao/jcvi/blob/main/src/jcvi/compara/synteny.py).

## minimap2 PAF

```bash
syntangle import-paf \
  --paf query_to_target.paf \
  --query-sizes query.fa.fai --target-sizes target.fa.fai \
  --query-species query --target-species target \
  --output imported_paf
```

One PAF alignment becomes one ribbon spanning its query and target intervals. The PAF strand controls target endpoint direction; coordinates remain zero-based boundaries. Sequence lengths must agree with the supplied sizes. CIGAR strings and within-alignment indels are not expanded, and separate alignment records are not merged into larger synteny blocks.

By default **all valid records are retained**, including secondary alignments and MAPQ 255 (unknown). To filter explicitly, add `--min-mapq 30` and/or `--min-alignment-length 10000`. A positive MAPQ threshold excludes unknown MAPQ. Alignment length means PAF column 11, including gaps, not just query span. Every excluded record's source line and reason is written to `excluded.tsv`; counts and thresholds are recorded in the audit. There is no hidden repeat filtering or best-hit selection.

Format reference: [minimap2 PAF specification](https://github.com/lh3/minimap2#output-format).

## Optimize the converted evidence

```bash
syntangle run --chromosomes imported_jcvi/chromosomes.tsv \
  --links imported_jcvi/links.tsv --output results/jcvi_full --figure
```

Or continue in R:

```r
library(SynTangleR)
x <- syntangle_read_tables("imported_jcvi")
fit <- syntangle_run(x, output = "results/jcvi_full", seconds = 600)
fit$stats
fit$pdf
```

The same R call accepts `imported_paf`. Table inputs already encode their coordinate system and initial signs, so omit `useOrder` and `invertTheseChrs`. Keep `import_audit.json` with downstream results; the R workflow copies it into the output directory.

For multi-species comparisons, combine pairwise converted link tables only after assigning globally unique link IDs and one complete chromosome table in the intended species row order. Links must join adjacent displayed species. Automatic multi-file merging is not part of these experimental importers.
