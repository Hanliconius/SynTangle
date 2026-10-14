function (refGenome, bed, clens, syntenyWeight) 
{
    ct_clust <- function(chr, ord, nclus) {
        if (any(duplicated(chr))) 
            stop("can't have duplicated chr values")
        if (length(ord) == 1) {
            names(ord) <- chr[1]
            return(ord)
        }
        else {
            if (length(ord) == 0) {
                return(1)
            }
            else {
                y <- data.matrix(ord)
                rownames(y) <- chr
                di <- dist(y)
                hc <- hclust(di, method = "ward.D2")
                nclus <- min(c(length(ord), nclus))
                ct <- cutree(hc, k = nclus)
                return(ct[chr])
            }
        }
    }
    ordByFun <- genome <- noAnchor <- refchr <- reford <- ord <- median <- plotOrd <- meanPos <- isArrayRep <- refChrOrder <- pSyn <- pName <- synPos <- meanOrd <- NULL
    tb <- subset(bed, !noAnchor & isArrayRep)[, c("genome", "chr", 
        "ord", "og")]
    tc <- clens[, c("genome", "chr", "ordByFun")]
    bedi <- merge(subset(tb, !duplicated(tb)), subset(tc, !duplicated(tc)), 
        by = c("genome", "chr"), allow.cartesian = T)
    bedr <- with(subset(bedi, genome == refGenome), data.table(refgenome = refGenome, 
        refchr = chr, refChrOrder = ordByFun, reford = ord, og = og))
    beda <- subset(bedi, genome != refGenome)
    bedm <- merge(subset(bedr, !duplicated(bedr)), subset(beda, 
        !duplicated(beda)), by = "og", allow.cartesian = T)
    setkey(bedm, refChrOrder, reford, genome, ord)
    bedm[, `:=`(reford, (1:.N)/.N), by = "genome"]
    bedm[, `:=`(pName = ordByFun/max(ordByFun), pSyn = reford/max(reford)), 
        by = "genome"]
    chro <- bedm[, list(meanPos = median(pSyn, na.rm = T), meanOrd = median(pName, 
        na.rm = T)), by = c("genome", "chr")]
    if (syntenyWeight <= 0) {
        nclus <- 1
    }
    else {
        if (syntenyWeight >= 1) {
            nclus <- Inf
        }
        else {
            nclus <- uniqueN(bedm$refchr)
        }
    }
    chro[, `:=`(synClus = ct_clust(ord = meanPos, chr = chr, 
        nclus = nclus)), by = "genome"]
    chro[, `:=`(synPos = mean(meanPos, na.rm = T)), by = c("genome", 
        "synClus")]
    setkey(chro, genome, synPos, meanOrd, chr)
    outa <- merge(clens, chro, by = c("genome", "chr"))
    outr <- subset(clens, genome == refGenome)
    outr[, `:=`(plotOrd, ordByFun)]
    setkey(outa, synPos, ordByFun, genome)
    outa[, `:=`(plotOrd, 1:.N), by = "genome"]
    outa <- outa[, colnames(outr), with = F]
    return(rbind(outr, outa))
}
