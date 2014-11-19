##uvSignif.r
##2014-10-27 dmontaner@cipf.es.com

##' @name uvSignif
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords significant enriched terms
##' @seealso \code{\link{uvGsa}}, \code{\link{uvPat}},
##' \code{\link{propagateGO}}, \code{\link{pval2index}}, \code{\link{goLeaves}}
##'
##' @title Filter significant terms in the univariate gene set analysis.
##'
##' @description 
##' Filters the rows in the data.frame returned by \code{\link{uvGsa}}
##' so that only the
##' enriched blocks are kept. 
##' 
##' @details
##' Works as \code{\link{goLeaves}} but without removing redundant terms.
##' 
##' @param gsaout data.frame; output from uvGsa.
##' @param cutoff p-value cutoff for considering significant a Gene Set.
##' @param pvalue p-value column to be used.
##' Default is named "padj" as in uvGsa output.
##' @param statistic name of the column containing the log odds ratio
##' from the uvGsa analysis.
##' @param sort if TRUE the output data.frame is ordered according to
##' significance.
##' @param verbose verbose
##' 
##' @return The input data.frame but keeping just the
##' 'significant' functional blocks.
##'
##' @examples
##' \dontrun{
##' res <- uvGsa (rindex, annotList)
##' uvSignif (res)
##' }
##' 
##' @export

uvSignif <- function (gsaout, cutoff = 0.05, pvalue = "padj",
                      statistic = "lor", verbose = TRUE, sort = TRUE) {
    
    ## real gsaout
    if (is.data.frame (gsaout) | is.matrix (gsaout)) {
        
        pat <- uvPat (gsaout, cutoff = cutoff, pvalue = pvalue,
                      statistic = statistic)
        
        touse <- pat %in% c (-1, 1)
        gsaout <- gsaout[touse, , drop = FALSE]
        
        ##sort
        if (sort) {
            myindex <- pval2index (pval = gsaout[,pvalue],
                                   sign = gsaout[,statistic],
                                   log = TRUE, verbose = verbose)
            orden <- order (myindex, decreasing = TRUE)
            gsaout <- gsaout[orden,]
        }
    } else {
        message ("gsaout must be a data.frame created by uvGsa")
    }
    ## RETURN
    gsaout
}
