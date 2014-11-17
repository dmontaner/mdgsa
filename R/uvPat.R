##uvPat.r
##2010-02-09 dmontaner@cipf.es
##2013-03-25 dmontaner@cipf.es


##' @name uvPat
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords univariate GSA pattern
##' @seealso \code{\link{uvGsa}}, \code{mdPat}
##'
##' @title Uni-Variate Gene Set Analysis Pattern Classification.
##' 
##' @description
##' Classifies significant patterns form a Uni-Variate Gene Set Analysis.
##' 
##' @details
##' Sign of the 'lor' and p-value are used to define functional blocks as
##' up-regulated, down-regulated or not enriched.
##' 
##' @param gsaout data.frame; output from uvGsa.
##' @param cutoff p-value cutoff for considering significant a Gene Set.
##' @param pvalue p-value column to be used. Default is named "padj" as in uvGsa output.
##' @param statistic name of the column containing the log odds ratio from the uvGsa analysis.
##' 
##' @return A numeric vector (values: -1, 0, 1) indicating relationship between the Gene Set and the ranking variable:
##' \describe{
##'   \item{1:}{indicates that the gene set is significantly associated to high values of the ranking statistic.}
##'   \item{-1:}{indicates that the gene set is significantly associated to low values of the ranking statistic.}
##'   \item{0:}{indicates that the gene set not related to the ranking statistic (no enrichment).}
##' }
##' 
##' @export
uvPat <- function (gsaout, cutoff = 0.05, pvalue = "padj", statistic = "lor") {
    res <- sign (gsaout[,statistic]) * as.numeric (gsaout[,pvalue] < cutoff)
    names (res) <- rownames (gsaout)
    res
}
