# complexion

[![Build Status](https://travis-ci.org/biobenkj/complexion.png?branch=master)](https://travis-ci.org/biobenkj/complexion)  [![codecov](https://codecov.io/gh/biobenkj/complexion/branch/master/graph/badge.svg)](https://codecov.io/gh/biobenkj/complexion)


Complexion is an R package for computing library complexity measures using rational function approximation as implemented in preseqR. Currently, it only supports reading in and processing BAM files. Though there are stubs that will eventually support splice junction bed files to compute complexity from.

### Quickstart

````
# list the BAMs
bams <- list.files(pattern="*.bam$")

# read in the BAMs aligned to hg38
# these were done with Ensembl annotations so instead of 'chr1' they are '1'
galignments.list <- lapply(bams, readBAM, genome="hg38", style="NCBI")

# estimate library complexity with 95% confidence intervals out to 100x the current library size
complexity.ests <- lapply(galignments.list, getEsts, withCI = TRUE)

# plot library complexity estimates with 95% confidence intervals
# note that the dotted line is the *current* complexity
plotComplexity(ests = complexity.ests[[1]], withCI = TRUE)
````