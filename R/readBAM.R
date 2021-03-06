#' pull in (filter and assemble as an object) some or all of the "useful" reads
#' in a paired-end sequencing run
#'
#' @param bam       character string, the BAM file to parse
#' @param genome    optional character string, the genome assembly to target
#' @param is.single whethr the input BAM is single-end or not
#' @param bamParams optional any parameters to pass through to ScanBamParam
#' @param which     optional a GRanges object with specific regions to extract
#' @param style     what style the assembly annotations are in (e.g. 'chr1' == UCSC or '1' == NCBI)
#' @param strand    what orientation is the first read in the pair (e.g. 'reverse')
#'
#' @return  a GenomicAlignmentPairs object
#'
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomicAlignments
#'
#' @export
#'

readBAM <- function(bam, genome=c("hg19", "hg38", "mm10", "mm9"),
                    is.single=FALSE, bamParams=NULL, which=NULL,
                    style="UCSC", strand=c("unstranded", "forward", "reverse"), ...) {
  
  getStdChromGRanges <- function(x) {
    gr <- keepStandardChromosomes(x, pruning.mode = "coarse")
    gr.prune <- keepSeqlevels(gr,
                              value = seqlevels(gr)[1:(which(seqlevels(gr) == "chrM") - 1)],
                              pruning.mode = "coarse")
    return(gr.prune)
  }
  
  if(is.null(bamParams)) {
    if (is.null(which)) {
      genome <- match.arg(genome)
      genome.gr <- switch(genome,
                          hg19 = data("hg19.gr", package="complexion"),
                          hg38 = data("hg38.gr", package="complexion"),
                          mm9 = data("mm9.gr", package="complexion"),
                          mm10 = data("mm10.gr", package="complexion"))
      which <- switch(genome,
                      hg19 = getStdChromGRanges(hg19.gr),
                      hg38 = getStdChromGRanges(hg38.gr),
                      mm9 = getStdChromGRanges(mm9.gr),
                      mm10 = getStdChromGRanges(mm10.gr))
    }
    #just in case we don't have UCSC convention
    if (style == "NCBI") seqlevelsStyle(which) <- "NCBI"
    if (style == "Ensembl") seqlevelsStyle(which) <- "Ensembl" 
    
    #make sure the index file is there
    if (!file.exists(paste0(bam, ".bai"))) {
      message("BAM index file doesn't exist for ", bam, ".")
      stop("Need to index the BAM before importing.")
    }

    #now check if we are in single-end mode
    if (is.single) {
      bamParams <- singleEndFilters(which=which, ...)
    } else {
      bamParams <- properPairedEndFilters(which=which, ...) 
      }
  }
  
  #set the strand
  strand <- match.arg(strand)
  strand <- switch(strand,
                   "unstranded"=0,
                   "stranded"=1,
                   "forward"=1,
                   "reverse"=2)
  
  if (is.single) {
    ga <- readGAlignments(bam, param=bamParams)
  } else {
    ga <- readGAlignmentPairs(bam, param=bamParams)
  }
  
  #set strandMode
  strandMode(ga) <- strand
  return(ga)
}

## not exported; used for preseq estimation
properPairedEndFilters <- function(which, ...) {
  
  ScanBamParam(what=c("rname","strand","pos","isize","mapq","qual","cigar"),
               flag=scanBamFlag(isProperPair=TRUE,
                                isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

singleEndFilters <- function(which, ...) {
  
  ScanBamParam(what=c("rname","strand","pos","mapq","qual","cigar"),
               flag=scanBamFlag(isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

## not exported; used for BAM filtering
uniquePairedEndFilters <- function(which, ...) {
  
  ScanBamParam(what=c("rname","strand","pos","isize","mapq","qual","cigar"),
               flag=scanBamFlag(isProperPair=TRUE,
                                isDuplicate=FALSE,
                                isNotPrimaryRead=FALSE,
                                isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

## for recovering spike-ins
spikeInFilters <- function(...) {
  
  ## grab the unmapped sequences for brute-force alignment to phiX spike-ins
  ScanBamParam(what=c("seq"), 
               flag=scanBamFlag(isUnmappedQuery=TRUE, 
                                isNotPassingQualityControls=FALSE, 
                                ...))
}
