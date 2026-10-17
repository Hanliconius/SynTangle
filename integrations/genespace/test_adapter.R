# Run from repository root. Python round-trip enabled by environment variable.
if (Sys.getenv('SYNTANGLE_TEST_R_INSTALLED')=='1') library(SynTangleR) else source('r/SynTangleR/R/flip_assist.R')
chr <- data.frame(genome=c('A','B'),chr=c('a','b'),length=100,x1=0,y1=c(1.1,2.1),y2=c(.9,1.9))
blk <- data.frame(genome1='A',chr1='a',genome2='B',chr2='b',
    startBp1=c(10,60),endBp1=c(20,70),startBp2=c(70,20),endBp2=c(60,10))
for (n in c('start','end')) for (i in 1:2) blk[[paste0(n,'Ord',i)]] <- blk[[paste0(n,'Bp',i)]]
rp <- list(plotData=list(sourceData=list(chromosomes=chr,blocks=blk)))
original <- serialize(rp,NULL)
root <- tempfile('syntangle-r-test-');dir.create(root)
export <- syntangle_export_genespace(rp,file.path(root,'export'),FALSE)
stopifnot(nrow(export$links)==2,export$links$start2[1]==70,export$links$end2[1]==60)
rankExport <- syntangle_export_genespace(rp,file.path(root,'rank'),TRUE)
stopifnot(identical(export$links,rankExport$links))
inv <- data.frame(genome='B',chr='b')
flipped <- syntangle_export_genespace(rp,file.path(root,'initial'),FALSE,inv)
stopifnot(identical(flipped$chromosomes$orientation,c(1L,-1L)))
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
bad <- rp;bad$plotData$sourceData$chromosomes <- chr[1,,drop=FALSE]
fail(syntangle_export_genespace(bad,file.path(root,'bad'),FALSE))
bad <- rp;bad$plotData$sourceData$blocks$startBp1[1] <- 20
fail(syntangle_export_genespace(bad,file.path(root,'zero'),FALSE))
# Exercise data.table-like subclasses without depending on data.table.
subclass <- rp
class(subclass$plotData$sourceData$blocks) <- c('data.table','data.frame')
class(subclass$plotData$sourceData$chromosomes) <- c('data.table','data.frame')
syntangle_export_genespace(subclass,file.path(root,'subclass'),FALSE)
exe <- Sys.getenv('SYNTANGLE_TEST_EXECUTABLE')
if (nzchar(exe)) {
    a <- syntangle_flip_assist(rp,FALSE,file.path(root,'flips with spaces'),executable=exe,figure=FALSE)
    stopifnot(a$stats$input_crossings==1,a$stats$output_crossings==0,a$stats$links==2,
              nrow(a$changedChromosomes)==1,nrow(a$layout)==2)
    b <- syntangle_flip_assist(rp,FALSE,file.path(root,'already'),invertTheseChrs=inv,executable=exe,figure=FALSE)
    stopifnot(b$stats$input_crossings==0,b$stats$output_crossings==0,nrow(b$changedChromosomes)==0,
              identical(b$invertTheseChrs,inv))
    c <- syntangle_run(rp,FALSE,file.path(root,'full'),seconds=3,executable=exe,figure=TRUE)
    stopifnot(c$stats$output_crossings==0,c$stats$links==2,file.exists(c$pdf),nrow(c$layout)==2)
    tableInput <- syntangle_read_tables(file.path(root,'export'))
    d <- syntangle_flip_assist(tableInput,output=file.path(root,'from tables'),executable=exe,figure=FALSE)
    stopifnot(d$stats$output_crossings==0,d$stats$links==2)
    # Shipped real-data tables exercise the installed R/Python path too.
    real <- syntangle_read_tables('examples/rhynchospora')
    for (workflow in c('flips','full')) {
        fun <- if (workflow=='flips') syntangle_flip_assist else syntangle_run
        result <- fun(real,output=file.path(root,paste0('rhynchospora_',workflow)),
                      executable=exe,seconds=20,figure=TRUE)
        stopifnot(result$stats$input_crossings==13,result$stats$output_crossings==9,
                  result$stats$links==37,nrow(result$layout)==10,file.exists(result$pdf))
    }
    print(a);print(c)
}
stopifnot(identical(serialize(rp,NULL),original))
unlink(root,recursive=TRUE)
cat('PASS: R schema, directed endpoints, original inversions and Python workflows\n')
