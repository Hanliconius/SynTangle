# Export the exact displayed chromosomes and directed region endpoints.
# Source this file; GENESPACE itself is not modified.
syntangle_export_genespace <- function(riparian, output, useOrder,
                                      invertTheseChrs = NULL) {
    if (length(useOrder) != 1L || is.na(useOrder) || !is.logical(useOrder))
        stop('useOrder must match the original plot_riparian call (TRUE or FALSE)')
    source <- if (!is.null(riparian$plotData)) riparian$plotData$sourceData else riparian$sourceData
    if (is.null(source)) stop('Expected plot_riparian result with plotData$sourceData')
    chromosomes <- as.data.frame(source$chromosomes)
    blocks <- as.data.frame(source$blocks)
    need <- function(x, columns) {
        if (!is.data.frame(x) || !all(columns %in% names(x)))
            stop('Missing source columns: ', paste(setdiff(columns,names(x)),collapse=', '))
    }
    need(chromosomes,c('genome','chr','length','x1','y1','y2'))
    suffix <- if (useOrder) 'Ord' else 'Bp'
    endpointNames <- paste0(rep(c('start','end'),2), suffix, rep(1:2,each=2))
    need(blocks,c('genome1','chr1','genome2','chr2',endpointNames))
    if (!nrow(chromosomes) || !nrow(blocks)) stop('No displayed chromosomes or links')
    key <- function(g,c) paste(g,c,sep='\034')
    keys <- key(chromosomes$genome,chromosomes$chr)
    if (anyDuplicated(keys)) stop('Duplicate displayed chromosome')
    numericCols <- c('length','x1','y1','y2')
    if (any(!is.finite(as.matrix(chromosomes[numericCols]))) || any(chromosomes$length<=0))
        stop('Invalid displayed lengths or coordinates')
    chromosomes <- chromosomes[order((chromosomes$y1+chromosomes$y2)/2,chromosomes$x1),]
    species <- unique(as.character(chromosomes$genome))
    # A species must occupy one horizontal row.
    if (any(vapply(split((chromosomes$y1+chromosomes$y2)/2,chromosomes$genome),function(x) length(unique(x))!=1L,logical(1))))
        stop('A species occupies multiple display rows')
    chr <- data.frame(species=as.character(chromosomes$genome),
        chromosome=as.character(chromosomes$chr),length=chromosomes$length,
        rank=ave(seq_len(nrow(chromosomes)),chromosomes$genome,FUN=seq_along),orientation=1L)
    keys <- key(chr$species,chr$chromosome)
    if (!is.null(invertTheseChrs) && nrow(invertTheseChrs)) {
        need(invertTheseChrs,c('genome','chr'))
        flips <- key(invertTheseChrs$genome,invertTheseChrs$chr)
        if (anyDuplicated(flips) || any(!flips %in% keys)) stop('Invalid original inversion list')
        chr$orientation[keys %in% flips] <- -1L
    }
    links <- data.frame(link_id=paste0('region_',seq_len(nrow(blocks))),
        species1=as.character(blocks$genome1),chromosome1=as.character(blocks$chr1),
        start1=blocks[[endpointNames[1]]],end1=blocks[[endpointNames[2]]],
        species2=as.character(blocks$genome2),chromosome2=as.character(blocks$chr2),
        start2=blocks[[endpointNames[3]]],end2=blocks[[endpointNames[4]]])
    for (side in 1:2) {
        idx <- match(key(links[[paste0('species',side)]],links[[paste0('chromosome',side)]]),keys)
        if (anyNA(idx)) stop('Source blocks reference chromosomes missing from the display; no links were discarded')
        a <- links[[paste0('start',side)]]; b <- links[[paste0('end',side)]]
        if (any(!is.finite(a)) || any(!is.finite(b)) || any(a==b) ||
            any(pmin(a,b)<0) || any(pmax(a,b)>chr$length[idx]))
            stop('Invalid directed endpoints; check useOrder and display extents')
    }
    if (any(abs(match(links$species1,species)-match(links$species2,species))!=1))
        stop('Links must join adjacent species rows')
    if (dir.exists(output) && length(list.files(output,all.files=TRUE,no..=TRUE)))
        stop('Output directory must be new or empty')
    dir.create(output,recursive=TRUE,showWarnings=FALSE)
    write.table(chr,file.path(output,'chromosomes.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    write.table(links,file.path(output,'links.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    writeLines(c(paste0('coordinate_mode=',suffix),
        'lengths=GENESPACE displayed extents; not necessarily full assembly lengths',
        paste0('GENESPACE_version=',if (requireNamespace('GENESPACE',quietly=TRUE)) as.character(packageVersion('GENESPACE')) else 'not installed')),file.path(output,'source.txt'))
    invisible(list(chromosomes=chr,links=links))
}

# Run the shared Python engine from an existing R session.
.syntangle_run_engine <- function(riparian, useOrder, invertTheseChrs,
                                 seconds, output, executable, figure, mode, extra=character()) {
    if (length(seconds)!=1L || !is.finite(seconds) || seconds<=0)
        stop('seconds must be finite and positive')
    if (length(executable)!=1L || !nzchar(executable)) stop('Provide the syntangle executable path')
    if (!nzchar(Sys.which(executable)) && !file.exists(executable))
        stop('Cannot find syntangle. Install the Python package or supply executable="/path/to/venv/bin/syntangle"')
    if (inherits(riparian,'syntangle_tables')) {
        if (!is.null(useOrder) || !is.null(invertTheseChrs))
            stop('Table inputs already encode coordinates and starting orientations; omit useOrder and invertTheseChrs')
        if (dir.exists(output) && length(list.files(output,all.files=TRUE,no..=TRUE)))
            stop('Output must be new or empty')
        dir.create(output,recursive=TRUE,showWarnings=FALSE)
        for (name in c('chromosomes.tsv','links.tsv','import_audit.json','excluded.tsv')) {
            src <- file.path(riparian$directory,name)
            if (file.exists(src) && !file.copy(src,file.path(output,name))) stop('Could not copy input ',name)
        }
        tables <- list(chromosomes=read.delim(file.path(output,'chromosomes.tsv'),colClasses='character'))
    } else {
        tables <- syntangle_export_genespace(riparian,output,useOrder,invertTheseChrs)
    }
    resultDir <- file.path(output,'result')
    args <- c(mode,'--chromosomes',shQuote(file.path(output,'chromosomes.tsv')),
        '--links',shQuote(file.path(output,'links.tsv')),'--output',shQuote(resultDir),
        '--seconds',as.character(seconds),if (figure) '--figure',extra)
    # R prepends its own shared-library directories to LD_LIBRARY_PATH.
    # Prefer the selected Python interpreter's lib directory (including venv
    # symlink resolution), so an unrelated system libpython cannot replace it.
    childEnv <- character()
    if (.Platform$OS.type=='unix') {
        launcher <- if (file.exists(executable)) executable else Sys.which(executable)
        launcher <- normalizePath(launcher)
        first <- readLines(launcher,n=1L,warn=FALSE)
        python <- if (length(first) && grepl('^#!.*[/]python[0-9.]*$',first)) sub('^#!','',first) else ''
        prefix <- if (nzchar(python) && file.exists(python)) dirname(dirname(normalizePath(python))) else dirname(dirname(launcher))
        lib <- file.path(prefix,'lib')
        if (dir.exists(lib)) childEnv <- paste0('LD_LIBRARY_PATH=',shQuote(paste(lib,Sys.getenv('LD_LIBRARY_PATH'),sep=':')))
    }
    status <- system2(executable,args,env=childEnv)
    if (status!=0) stop('SynTangle failed with exit status ',status,'; exported inputs retained in ',output)
    layout <- read.delim(file.path(resultDir,'layout.tsv'),colClasses='character',check.names=FALSE)
    negative <- layout$orientation=='-1'
    result <- list(invertTheseChrs=data.frame(genome=layout$species[negative],chr=layout$chromosome[negative]),
        changedChromosomes=read.delim(file.path(resultDir,'flips.tsv'),colClasses='character'),
        layout=layout,inputLayout=tables$chromosomes,
        stats=read.delim(file.path(resultDir,'summary.tsv'),check.names=FALSE),
        pairStats=read.delim(file.path(resultDir,'pair_stats.tsv'),check.names=FALSE),
        audit=file.path(resultDir,'stats.json'),
        pdf=if (figure) file.path(resultDir,paste0(gsub('-','_',mode),'_comparison.pdf')) else NULL,
        mode=mode)
    class(result) <- 'syntangle_result'
    saveRDS(result,file.path(output,'syntangle_result.rds'))
    result
}

syntangle_flip_assist <- function(riparian, useOrder=NULL, output, invertTheseChrs=NULL,
                                 seconds=60, executable='syntangle', figure=TRUE) {
    .syntangle_run_engine(riparian,useOrder,invertTheseChrs,seconds,output,executable,figure,'flip-assist')
}

syntangle_run <- function(riparian, useOrder=NULL, output, invertTheseChrs=NULL,
                         seconds=600, executable='syntangle', figure=TRUE,
                         ribbonObjective=c('none','crossing-cap','combined'),
                         twistMetric=c('width','count'), weight=1, ribbonSeconds=600, seed=1L) {
    ribbonObjective <- match.arg(ribbonObjective); twistMetric <- match.arg(twistMetric)
    if (length(weight)!=1L || !is.finite(weight) || weight<0) stop('Invalid weight')
    if (length(ribbonSeconds)!=1L || !is.finite(ribbonSeconds) || ribbonSeconds<=0) stop('Invalid ribbonSeconds')
    if (length(seed)!=1L || !is.finite(seed) || seed!=as.integer(seed)) stop('Invalid seed')
    extra <- c('--ribbon-objective',ribbonObjective,'--twist-metric',twistMetric,
        '--weight',as.character(weight),'--ribbon-seconds',as.character(ribbonSeconds),'--seed',as.character(seed))
    .syntangle_run_engine(riparian,useOrder,invertTheseChrs,seconds,output,executable,figure,'run',extra)
}

print.syntangle_result <- function(x,...) {
    cat('SynTangle',x$mode,':',x$stats$input_crossings,'->',x$stats$output_crossings,'midpoint crossings\n')
    cat(nrow(x$changedChromosomes),'whole-chromosome flips;',x$stats$links,'links retained\n')
    cat('Status:',x$stats$status,'; independent certificate: no\n')
    invisible(x)
}

# Read the canonical tables, including outputs from experimental format importers.
syntangle_read_tables <- function(directory) {
    if (!all(file.exists(file.path(directory,c('chromosomes.tsv','links.tsv')))))
        stop('Directory must contain chromosomes.tsv and links.tsv')
    structure(list(directory=normalizePath(directory)),class='syntangle_tables')
}
