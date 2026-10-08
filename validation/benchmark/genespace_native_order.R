# Execute the installed GENESPACE ordering code, without its rendering stack.
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2L, 3L))) stop("Usage: genespace_native_order.R INPUT OUTPUT [--skip-incomplete-variants]")
skipIncomplete <- length(args) == 3L && args[[3]] == "--skip-incomplete-variants"
if (length(args) == 3L && !skipIncomplete) stop("Unknown optional argument")
library(data.table)
ns <- asNamespace("GENESPACE")
version <- as.character(packageVersion("GENESPACE"))
if (version != "1.3.1") stop("Pilot validated for GENESPACE 1.3.1; found ", version)

# This function is nested inside riparian_engine(), not exported. Extract the
# exact assignment from the installed function body; never reimplement it.
engine <- get("riparian_engine", envir = ns)
assignments <- Filter(function(x) {
  is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], as.name("pull_synChrOrd"))
}, as.list(body(engine))[-1L])
if (length(assignments) != 1L) stop("Installed ordering function shape changed")
scope <- new.env(parent = ns)
eval(assignments[[1]], envir = scope)
native_order <- scope$pull_synChrOrd

input <- args[[1]]
output <- args[[2]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
bed <- fread(file.path(input, "bed.tsv"))
clens <- fread(file.path(input, "clens.tsv"))
species <- fread(file.path(input, "species.tsv"))$species_id
required <- c("genome", "chr", "ord", "og", "noAnchor", "isArrayRep")
if (!all(required %in% names(bed))) stop("Invalid public-evidence bed")
if (anyDuplicated(bed[, .(genome, og)])) stop("Pilot requires unique homology per species")

rows <- list()
variantAudit <- list()
for (reference in species) {
  for (weight in c(1, 0.5)) {
    id <- sprintf("%s_w%s", reference, weight)
    started <- proc.time()[["elapsed"]]
    order <- tryCatch(native_order(reference, copy(bed), copy(clens), weight),
                      error = function(e) {
                        if (!skipIncomplete) stop(e)
                        structure(list(message = conditionMessage(e)), class = "native_order_error")
                      })
    seconds <- proc.time()[["elapsed"]] - started
    expected <- paste(clens$genome, clens$chr, sep = "::")
    failed <- inherits(order, "native_order_error")
    returned <- if (failed) character() else paste(order$genome, order$chr, sep = "::")
    missing <- setdiff(expected, returned)
    extra <- setdiff(returned, expected)
    complete <- !failed && length(returned) == length(expected) &&
      !anyDuplicated(returned) && length(missing) == 0L && length(extra) == 0L
    variantAudit[[length(variantAudit) + 1L]] <- data.table(
      variant = id, reference = reference, synteny_weight = weight,
      complete = complete, expected_chromosomes = length(expected),
      returned_chromosomes = length(returned), missing = paste(missing, collapse = ";"),
      extra = paste(extra, collapse = ";"),
      error = if (failed) order$message else "", ordering_seconds = seconds)
    fwrite(rbindlist(variantAudit), file.path(output, "native_variant_audit.tsv"), sep = "\t")
    if (!complete) {
      if (!skipIncomplete) stop("GENESPACE omitted chromosomes; evidence will not be silently dropped")
      message("SKIP ", id, ": incomplete native result; missing=", length(missing),
              "; no incomplete layout enters the comparison")
      next
    }
    if (anyDuplicated(order[, .(genome, chr)])) stop("Duplicate output chromosomes")
    order[, `:=`(variant = id, reference = reference, synteny_weight = weight,
                 ordering_seconds = seconds)]
    rows[[length(rows) + 1L]] <- order[, .(
      variant, reference, synteny_weight, genome, chr, plotOrd, ordering_seconds
    )]
  }
}
if (length(rows) == 0L) stop("No complete GENESPACE ordering variant; inspect native_variant_audit.tsv")
fwrite(rbindlist(rows), file.path(output, "native_orders.tsv"), sep = "\t")
writeLines(deparse(native_order), file.path(output, "installed_ordering_function.R"))
writeLines(c(paste("GENESPACE", version), R.version.string,
             paste("data.table", packageVersion("data.table"))),
           file.path(output, "versions.txt"))
