# Execute the installed GENESPACE ordering code, without its rendering stack.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: genespace_native_order.R INPUT OUTPUT")
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
for (reference in species) {
  for (weight in c(1, 0.5)) {
    started <- proc.time()[["elapsed"]]
    order <- native_order(reference, copy(bed), copy(clens), weight)
    seconds <- proc.time()[["elapsed"]] - started
    if (nrow(order) != nrow(clens) ||
        !setequal(paste(order$genome, order$chr), paste(clens$genome, clens$chr))) {
      stop("GENESPACE omitted chromosomes; evidence will not be silently dropped")
    }
    if (anyDuplicated(order[, .(genome, chr)])) stop("Duplicate output chromosomes")
    id <- sprintf("%s_w%s", reference, weight)
    order[, `:=`(variant = id, reference = reference, synteny_weight = weight,
                 ordering_seconds = seconds)]
    rows[[length(rows) + 1L]] <- order[, .(
      variant, reference, synteny_weight, genome, chr, plotOrd, ordering_seconds
    )]
  }
}
fwrite(rbindlist(rows), file.path(output, "native_orders.tsv"), sep = "\t")
writeLines(deparse(native_order), file.path(output, "installed_ordering_function.R"))
writeLines(c(paste("GENESPACE", version), R.version.string,
             paste("data.table", packageVersion("data.table"))),
           file.path(output, "versions.txt"))
