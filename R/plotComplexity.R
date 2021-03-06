#' plot results from preseq (complexity estimation around features)
#' TODO: use the ... to allow for comparison of multiple libraries
#'
#' @param   x       the thing whose complexity is to be plotted
#' @param   ...     other parameters to pass on to getEsts()
#' @param   ests    precomputed complexity estimates (speeds up multiple runs) 
#' @param   withCI  plot confidence intervals for the estimates?  (FALSE)
#' 
#' @return  invisibly, the complexity estimates corresponding to the plot
#' 
#' @seealso preseqR
#' 
#' @import ggplot2
#' 
#' @export
#'

plotComplexity <- function(x, ..., ests=NULL, withCI=FALSE) {
  if (!is.null(ests)) { 
    if (withCI & !all(c("lb","ub") %in% names(ests))) {
      stop("CIs requested, but provided estimates do not contain them!")
    }
    message("Using supplied estimates...")
    res <- ests
  } else if (is(x, "data.frame")) {
    res <- getEsts(x, ..., withCI=withCI)
  } else if (is(x, "GRanges") | is(x, "GAlignmentPairs")) {
    res <- getEsts(getComplexity(x), ..., withCI=withCI)
  } else {
    stop("Don't know what to do with a", class(x))
  }
  
  if (is.null(res)) {
    stop("Complexity estimate failed. You might try estimating from splice junctions instead.")
  } else {
    message("Plotting complexity curve", 
            ifelse(withCI, " with confidence intervals...", "..."))
    if (withCI) { 
      epsilon <- 1000000
      ggplot(res, aes(x = reads/epsilon, y = frags/epsilon)) +
        geom_line() +
        geom_point() +
        geom_ribbon(aes(ymin = lb/epsilon, ymax = ub/epsilon, alpha = 0.2)) +
        geom_vline(xintercept = res$reads[1]/epsilon, linetype = "dashed") + #Add a line for existing library complexity
        xlab("Read depth (M)") +
        ylab("Unique fragments (M)") +
        theme(
          legend.title = element_blank(),
          legend.position = "None"
        )
    } 
    else {
      epsilon <- 1000000
      ggplot(res, aes(x = reads/epsilon, y = frags/epsilon)) +
        geom_line() +
        geom_point() +
        geom_vline(xintercept = res$reads[1]/epsilon, linetype = "dashed") + #Add a line for existing library complexity
        xlab("Read depth (M)") +
        ylab("Unique fragments (M)") +
        theme(
          legend.title = element_blank()
        )
    }
  }
}