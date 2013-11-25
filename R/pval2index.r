##pval2index.r
##2013-05-01 dmontaner@cipf.es


##' @name pval2index
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords p-value ranking index
##' @seealso \code{\link{indexTransform}}
##'
##' @title Transform p-values in into a ranking index.
##' 
##' @description
##' After a genomic test, p-values are numerical indexes which account for certain biological characteristic.
##' By definition p-values are bounded between zero and one, but this may not be suitable as an index.
##' Moreover, p-values are always derived form a statistic which sign may be important.
##' The function helps transforming the p-value and its associated statistic into a ranking index.
##' 
##' @details
##' By default the transformation is (-1) * log (pval) * sign (sign).
##' When log = FALSE the transformation is (1 - pval) * sign (sign).
##'
##' @param pval a vector or matrix of p-values.
##' @param sign a vector or matrix of signs associated to the p-values.
##' @param log = TRUE
##' @param verbose verbose
##'
##' @return A transformed index. A vector or matrix, depending on the input parameters.
##'
##' @examples 
##' my.statistic <- rnorm (1000)
##' my.pvalue <- 2 * pnorm (my.statistic)
##' my.pvalue[my.pvalue > 1] <- 2 - my.pvalue[my.pvalue > 1]
##' 
##' index <- pval2index (pval = my.pvalue, sign = my.statistic)
##' 
##' #par (mfrow = c (1,2))
##' #plot (my.statistic, my.pvalue)
##' #plot (my.statistic, index)
##' 
##' @export
pval2index <- function (pval, sign, log = TRUE, verbose = TRUE) {
  
  pval.cero <- pval == 0
  if (sum (pval.cero) & verbose) {
    cat ("\n", "Some p-values are strictly zero; infinite values will be returned.", "\n", fill = TRUE)
  }
  
  sign.cero <- sign == 0
  if (sum (sign.cero) & verbose) {
    cat ("\n", "Some p-values are strictly zero; infinite values will be returned.", "\n", fill = TRUE)
  }
  
  if (log) {
    res <- (-1) * log (pval) * sign (sign)
  } else {
    res <- (1 - pval) * sign (sign)
  }
  
  return (res)
}
