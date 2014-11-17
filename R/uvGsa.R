##uvGsa.r
##2009-11-16 dmontaner@cipf.es
##2013-03-25 dmontaner@cipf.es

## To Do:
##  CHECK na.action
##  CHECK weights

##' @name uvGsa
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases univariateGsa
##' 
##' @keywords univariate GSA gene set
##' @seealso \code{\link{mdGsa}}, \code{\link{uvPat}}, \code{glm.fit}, \code{p.adjust}
##'
##' @title Uni-Variate Gene Set Analysis.
##' 
##' @description
##' Performs a Uni-Variate Gene Set Analysis using a logistic regression model.
##' 
##' @details
##' 'index' may also be a numerical \code{matrix} or \code{data.frame}.
##' If such a matrix has more than one column, the ranking index is taken form the first one.
##' The remaining columns are used as covariates to correct for within the analysis.
##' 
##' Default p-value correction is "BY".
##' 
##' @param index ranking index, generally a numerical named vector.
##' @param annot an annotation list.
##' @param p.adjust.method p-value adjustment method for multiple testing.
##' @param family see \code{glm.fit}.
##' @param fulltable if TRUE, 'sd', 't' and 'convergence' indicator from the glm fit are included in the output.
##' @param verbose verbose.
##' @param verbosity integer indicating which iterations should be indicated if verbose = TRUE.
##' @param \dots further arguments to be pasted to \code{glm.fit}, for instance 'weights'.
##' 
##' @return A \code{data.frame} with a row for each Gene Set or block.
##' Columns are:
##' \describe{
##'   \item{\code{N}:}{number of genes annotated to the Gene Set.}
##'   \item{\code{lor}:}{log Odds Ratio estimated for the Gene Set.}
##'   \item{\code{pval}:}{p-values associated to each log Odds Ratio.}
##'   \item{\code{padj}:}{adjusted p-values.}
##'   \item{\code{sd}:}{standard deviations associated to each log Odds Ratio.}
##'   \item{\code{t}:}{t statistic associated to each log Odds Ratio.}
##' }
##' 
##' @export
uvGsa <- function (index, annot, p.adjust.method = "BY", family = quasibinomial(),
                   verbose = TRUE, verbosity = 100, fulltable = FALSE, ...) {

  
  ## INPUT FORMATTING ##########################################################
  
  ##GENE IDS
  if (is.data.frame (index)) {
    index <- as.matrix (index)
  }
  if (is.matrix (index)) {
    genes <- rownames (index)
    rownames (index) <- NULL
  }
  if (is.vector (index)) {
    genes <- names (index)
    names (index) <- NULL  ##somehow quicker ???
  }
  if (is.null (genes)) {
    stop ("no geneIds found. Check names or rownames in 'index'")
  }
  
  ##BLOCK IDS
  blocks <- names (annot)
  if (is.null (blocks)) {
    stop ("unnamed 'annot'")
  }

  
  ## INTERNALS #################################################################
  
  ##results matrix
  res <- matrix (NA, nrow = length (blocks), ncol = 6)
  rownames (res) <- blocks
  colnames (res) <- c("N", "lor", "sd", "t", "pval", "conv") # "error"

  ##store time
  t0 <- proc.time ()
  
  ##set counter
  counter <- 0
  if (verbose) {
    message ("Analyzed blocks:")
  }

  
  ## ALGORITHM #################################################################
  
  X <- cbind (rep (1, times = length (genes)), index) ##index may be a matrix or a vector
  colnames (X) <- NULL
  for (bl in blocks) {
    B <- as.numeric (genes %in% annot[[bl]])
    res.glm <- glm.fit (x = X, y = B, family = family, ...)
    res.sum <- summary.glm (res.glm)
    res[bl,] <- c (sum (B), res.sum$coefficients[2,], res.glm$converged)
    ##
    if (verbose) {
      counter <- counter + 1
      if (counter %% verbosity == 0) {
        message (counter, ", ", appendLF = FALSE)
      }
      if (counter %% (10*verbosity) == 0) {
        message ("\n")
      }
    }
  }
  ## ###########################################################################

  
  ##Time
  t1 <- proc.time ()
  if (verbose) {
    message ("time in seconds:")
    print (t1-t0)
  }
  
  ##Convergence
  if (verbose) {
    if (any (res[,"conv"] == 0)) {
      warning (paste ("The analysis did not converge for some blocks.",
                      "You may re-run uvGsa using 'fulltable = TRUE' to find them.",
                      sep = "\n"))
    }
  }

  
  ## p-value ADJUSTMENT #######################################################
  
  res <- cbind (res, padj = p.adjust (res[,"pval"], method = p.adjust.method))

  
  ##OUTPUT ####################################################################

  ##reorder columns
  res <- res[,c("N", "lor", "pval", "padj", "sd", "t", "conv")] # "error"
  
  ##remove some columns
  if (!fulltable) {
    res <- res[,c("N", "lor", "pval", "padj")]
  }
  
  ##format data.frame
  res <- as.data.frame (res)

  ##return
  res
}


