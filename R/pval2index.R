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
##' After a genomic test, p-values are numerical indexes which account for
##' certain biological characteristic.
##' By definition p-values are bounded between zero and one,
##' but this may not be suitable as an index.
##' Moreover, p-values are always derived form a statistic which sign may be
##' important.
##' The function helps transforming the p-value and its associated statistic
##' into a ranking index.
##' 
##' @details
##' The default transformation is (-1) * log (pval) * sign (sign).
##' When log = FALSE the transformation is (1 - pval) * sign (sign).
##'
##' If \code{sign} is missing all p-values are associated wit a positive sign.
##'
##' Missing values are allowed and return NA values.
##' 
##' An \code{offset} may be provided to replace p-values equal to zero
##' when \code{log = TRUE}.
##' In such way \strong{infinite} values are not generated.
##' If the \code{offset} parameter is not provided,
##' the minimum p-value other than zero is used for the replacement.
##' You can explicitly specify \code{offset = 0} if you want \code{Inf}
##' values to be returned.
##'
##' By default the names of the output vector (or row names in a matrix)
##' are those of \code{pval} or \code{sign}.
##' If \code{names} is provided, then it is used instead.
##'
##' @param pval a vector or matrix of p-values.
##' @param sign a vector or matrix of signs associated to the p-values.
##' @param names a character vector of the names of the features.
##' @param log = TRUE
##' @param offset value used to replace p-values equal to zero
##' @param verbose verbose
##'
##' @return A transformed index.
##' A vector or matrix, depending on the input parameters.
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
##'
##' ## Zero p-values
##' p <- c (0:10)/10
##' p
##' pval2index (p)
##' pval2index (p, offset = 0)
##' pval2index (p, offset = 0.000001)
##'
##' ## Missing p-values
##' p <- c(0:10, NA)/10
##' p
##' pval2index (p)
##' pval2index (p, offset = 0)
##' pval2index (p, offset = 0.000001)
##' pval2index (p, log = FALSE)
##' pval2index (p, offset = 0, log = FALSE)
##'
##' ## Matrix
##' p <- matrix (c(0:10, NA)/10, ncol = 3)
##' p
##' pval2index (p)
##' pval2index (p, offset = 0)
##' pval2index (p, offset = 0.000001)
##' pval2index (p, log = FALSE)
##' pval2index (p, offset = 0, log = FALSE)
##' 
##' @export
pval2index <- function (pval, sign, names = NULL, log = TRUE,
                        offset, verbose = TRUE) {
    
    if (missing (pval)) {
        stop ("pval is missing with no default")
    }
    
    if (missing (sign)) {
        message ("sign is missing. All signs will be considered as positive")
        sign <- rep (1, times = length (pval))
    }
    
    ## zero p-values and log transformation
    pval.cero <- pval == 0
    if (any (pval.cero, na.rm = TRUE) & log) {
        if (missing (offset)) {
            touse <- !(pval.cero | is.na (pval)) #note: (NA | TRUE) is TRUE
            offset <- min (pval[touse])
        }
        pval[which (pval.cero)] <- offset
    }
    
    sign.cero <- sign == 0
    if (any (sign.cero, na.rm = TRUE) & verbose) {
        message ("Some sign statistics are zero; zero values will be returned.")
    }
    
    if (log) {
        res <- (-1) * log (pval) * sign (sign)
    } else {
        res <- (1 - pval) * sign (sign)
    }
    
    if (!is.null (names)) {
        names (res) <- names
    }

    ## OUTPUT
    res
}
