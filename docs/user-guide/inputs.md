# Input tables

Use tab-separated UTF-8 files with headers. Empty tables, missing values and ambiguous links are rejected. Species row order is the order of first appearance in `chromosomes.tsv`; it remains fixed during optimization.

## Chromosomes

Required columns: `species`, `chromosome`, `length`, `rank`, `orientation`.

- `(species, chromosome)` uniquely identifies a chromosome. IDs must not contain `:`.
- `length` is finite and positive, in the same units as link endpoints.
- `rank` is a positive integer, unique within each species, specifying left-to-right order.
- `orientation` is `1` or `-1`. Coordinates always refer to the original chromosome, irrespective of displayed orientation.
- Retain chromosomes without links if they are part of the intended display.

## Directed links

Required columns: `link_id`, `species1`, `chromosome1`, `start1`, `end1`, `species2`, `chromosome2`, `start2`, `end2`.

Each row is one explicit ribbon between adjacent displayed species. IDs must be unique. Both endpoints must reference chromosomes in the chromosome table. All endpoint values must be finite and within `[0,length]`; each interval must have nonzero width.

`start1` corresponds to `start2`; `end1` corresponds to `end2`. Reverse the endpoint direction on one side to encode an inverted ribbon. **Do not sort both sides independently**: doing so erases relative orientation. A whole-chromosome flip transforms both endpoints by `x -> length - x`.

These are continuous boundary coordinates. For a conventional 1-based inclusive GFF interval `[s,e]`, use boundaries `[s-1,e]` before choosing its direction. GENESPACE exports use its native displayed Bp or Ord coordinates and displayed extents; do not independently substitute assembly lengths.

The bundled [example](../../examples/quickstart) illustrates all required columns. Original tables are hashed, and the normalized fixed evidence is saved as `fixture.json`.

## JSON fixtures and saved layouts

Advanced users may supply `--fixture input.json` instead of the two tables. See [the data model](../data_model.md). Public workflows require known block strands and unambiguous homology occurrences within species. Use explicit pairwise link tables for repeated copies.

`--layout layout.tsv` overrides the starting display and requires a complete table of `species`, `chromosome`, `rank`, `orientation`. It must preserve all chromosome identities and satisfy supplied orientation constraints.
