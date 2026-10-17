# Troubleshooting

**R cannot find `syntangle`.** Install the Python package into a virtual environment and pass the absolute executable path, for example `executable="/home/me/SynTangle/.venv/bin/syntangle"`. RStudio need not inherit the environment activated in another terminal.

**Python works in the terminal but fails inside R.** Supply the absolute virtual-environment executable path. The R launcher prioritizes that interpreter's shared-library directory to avoid loading an unrelated system Python runtime through R's library search path. Include both Python and R versions if the problem persists.

**Missing source-data columns.** Confirm the input is the return value of `plot_riparian()`, not only its ggplot object or `gsParam`. Report the GENESPACE version and `names(rp$plotData$sourceData$blocks)`.

**Missing displayed chromosomes.** The source blocks reference chromosomes absent from the GENESPACE display. Correct the original display and its ordering/filtering arguments. The public adapter rejects this mismatch rather than inventing an order or dropping links.

**Endpoint outside chromosome extent.** Check that `useOrder` matches the original plot. Do not mix gene ranks with base pairs, or replace displayed extents with unrelated assembly lengths.

**Unexpected starting orientation.** Pass the original `invertTheseChrs` list to the R helper. When replotting, replace it with the returned list; do not concatenate lists.

**Twisted ribbons remain after fewer crossings.** Midpoint crossings and ribbon twists are different quantities. Use optional crossing-cap refinement to improve twists without increasing the crossing score reached by the first stage.

**No improvement or nonzero gap.** This can reflect a time limit, a local optimum, or an already optimal layout. A nonzero lower-bound gap is unresolved; it is not proof that the layout is suboptimal.

**Output directory exists.** Choose a new directory for each run. This protects earlier statistics and figures.

For other failures, open an issue with the command or R call, versions, error text and a small shareable example. A successful file collector or ZIP creation is not evidence that all analyses inside it succeeded.
