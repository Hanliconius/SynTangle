args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop('Usage: rhynchospora_native_regions.R RUN MCScanX_DIRECTORY')
run <- normalizePath(args[1]); mcscan <- normalizePath(args[2])
suppressPackageStartupMessages(library(GENESPACE))
suppressPackageStartupMessages(library(data.table))
if(as.character(packageVersion('GENESPACE')) != '1.3.1') stop('Requires GENESPACE 1.3.1')
gids <- c('R_tenuis','R_austrobrasiliensis','R_breviuscula')
gpar <- init_genespace(wd=file.path(run,'genespace'),genomeIDs=gids,ploidy=rep(1,3),
 path2orthofinder='orthofinder',path2diamond='diamond',path2mcscanx=mcscan,
 nCores=as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','4')))
cat('GENESPACE DISCOVERY START\n');flush.console()
gpar <- run_genespace(gpar,makePairwiseFiles=TRUE)
saveRDS(gpar,file.path(run,'genespace_parameters.rds'))
for(regions in c(TRUE,FALSE)) {
 label <- if(regions) 'regions' else 'blocks'
 cat(sprintf('NATIVE PLOT %s START\n',label));flush.console()
 out <- plot_riparian(gsParam=gpar,refGenome='R_breviuscula',genomeIDs=gids,
  palette=colorRampPalette(c('#F9AC60','#307BB5','#D61F27','#ADD8E7','#FCF7BF')),
  braidAlpha=.75,chrExpand=.75,chrLabFontSize=8,chrBorderLwd=.3,chrBorderCol='black',
  useOrder=FALSE,useRegions=regions,forceRecalcBlocks=TRUE,minChrLen2plot=0,
  invertTheseChrs=data.frame(genome='R_breviuscula',chr='Chr3_h1'),
  customRefChrOrder=c('Chr2_h1','Chr5_h1','Chr1_h1','Chr4_h1','Chr3_h1'),
  pdfFile=file.path(run,paste0('native_',label,'.pdf')),scalePlotHeight=2)
 saveRDS(out,file.path(run,paste0('native_',label,'.rds')))
 fwrite(out$plotData$sourceData$blocks,file.path(run,paste0('native_',label,'.tsv')),sep='\t')
 fwrite(out$plotData$sourceData$chromosomes,file.path(run,paste0('native_',label,'_chromosomes.tsv')),sep='\t')
}
writeLines('PASS',file.path(run,'NATIVE_COMPLETE'))
