\name{syntangle_flip_assist}
\alias{syntangle_flip_assist}
\alias{syntangle_run}
\alias{syntangle_export_genespace}
\title{Optimize an Existing GENESPACE Chromosome Display}
\description{Call the SynTangle Python engine from R on the source data returned by GENESPACE. Internal coordinates and all supplied links are preserved.}
\usage{
syntangle_flip_assist(riparian, useOrder, output, invertTheseChrs = NULL,
    seconds = 60, executable = "syntangle", figure = TRUE)
syntangle_run(riparian, useOrder, output, invertTheseChrs = NULL,
    seconds = 600, executable = "syntangle", figure = TRUE,
    ribbonObjective = c("none", "crossing-cap", "combined"),
    twistMetric = c("width", "count"), weight = 1,
    ribbonSeconds = 600, seed = 1L)
syntangle_export_genespace(riparian, output, useOrder,
    invertTheseChrs = NULL)
}
\arguments{
\item{riparian}{The returned object from GENESPACE's plot_riparian function.}
\item{useOrder}{Logical matching the original plot: TRUE for gene ranks, FALSE for base-pair coordinates. Required explicitly.}
\item{output}{New or empty output directory.}
\item{invertTheseChrs}{The original plot's inversion data frame with columns genome and chr, or NULL.}
\item{seconds}{Positive search budget in seconds.}
\item{executable}{SynTangle executable name or absolute path to its Python virtual environment executable.}
\item{figure}{Write a vector before/after PDF.}
\item{ribbonObjective}{Optional refinement: none, crossing-cap or combined.}
\item{twistMetric}{Width-weighted or unweighted twist objective.}
\item{weight}{Nonnegative twist weight for the combined objective.}
\item{ribbonSeconds}{Additional refinement budget in seconds.}
\item{seed}{Integer search seed.}
}
\value{Optimization functions return a syntangle_result list containing layout, inputLayout, stats, pairStats, changedChromosomes, invertTheseChrs, audit and pdf. The inversion list is absolute and replaces the original list. The export function invisibly returns the written chromosome and link tables.}
\details{Flip assistance preserves chromosome order and performs whole-chromosome flips only. Full optimization also permits chromosome permutations. The PDF from SynTangle renders the exact resulting layout. Passing only the inversion list back to GENESPACE does not apply a full optimized chromosome order. Existing inversions must be supplied because GENESPACE source blocks precede the display inversion transform. Numerical bound closure is distinct from independent certification.}
\examples{
\dontrun{
library(GENESPACE)
rp <- plot_riparian(gsParam = gsParam, useOrder = FALSE)
fit <- syntangle_flip_assist(rp, useOrder = FALSE, output = "flips")
fit$stats
full <- syntangle_run(rp, useOrder = FALSE, output = "full")
full$layout
}
}
