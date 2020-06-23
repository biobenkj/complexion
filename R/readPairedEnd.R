#' pull in (filter and assemble as an object) some or all of the "useful" reads
#' in a paired-end sequencing run
#'
#' @param bam       character string, the BAM file to parse
#' @param genome    optional character string, the genome assembly to target
#' @param bamParams optional any parameters to pass through to ScanBamParam
#' @param which     optional a GRanges object with specific regions to extract
#'
#' @return  a GenomicAlignmentPairs object
#'
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomicAlignments
#'
#' @export
#'

readPairedEnd <- function(bam, genome=c("hg19", "hg38", "mm10", "mm9"),
                          bamParams=NULL, which=NULL, ...) {
  
  getStdChromGRanges <- function(x) {
    ## ONLY works if chromosomes are properly ordered as in OrganismDbi
    as(seqinfo(x), "GRanges")[ 1:(which(seqlevels(x) == "chrM") - 1) ] 
  }
  
  if(is.null(bamParams)) {
    if (is.null(which)) {
      genome <- match.arg(genome)
      which <- switch(genome,
                      hg19 = getStdChromGRanges(data("hg19.gr", package="complexion")),
                      hg38 = getStdChromGRanges(data("hg38.gr", package="complexion")),
                      mm9 = getStdChromGRanges(data("mm9.gr", package="complexion")),
                      mm10 = getStdChromGRanges(data("mm10.gr", package="complexion")))
    } 
    bamParams <- properPairedEndFilters(which=which, ...)
  }
  
  readGAlignmentPairs(bam, param=bamParams)
  
}

## not exported; used for preseq estimation
properPairedEndFilters <- function(which, ...) {
  
  ScanBamParam(what=c("rname","strand","pos","isize","mapq"),
               flag=scanBamFlag(isProperPair=TRUE,
                                isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

## not exported; used for BAM filtering
uniquePairedEndFilters <- function(which, ...) {
  
  ScanBamParam(what=c("rname","strand","pos","isize","mapq"),
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