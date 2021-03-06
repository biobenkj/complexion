% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getEsts.R
\name{getEsts}
\alias{getEsts}
\title{Complexity estimation main function.}
\usage{
getEsts(xx, withCI = FALSE, ...)
}
\arguments{
\item{xx}{The fragments or sample of fragments}

\item{withCI}{Have preseq compute 95 percent confidence intervals for plots?}

\item{...}{Other arguments to pass on to preseqR}
}
\value{
\preformatted{   A data frame with results
}
}
\description{
Estimate the complexity of a library or sample based on unique fragments
using Daley and Smith's implementation of Good-Toulmin rational function
approximation to solve the missing species problem.
}
\details{
Original functions in preseqR v2.0.1.1 for this were:
preseqR.rfa.curve and preseqR.rfa.species.accum.curve

The new functions as of the version 4.0.0 are:

ds.rSAC == preseqR.rfa.curve
ds.rSAC.bootstrap == preseqR.rfa.species.accum.curve

The new functions return generators that can have data passed to them
instead of returning a data frame as in version 2.0.1.1.
}
