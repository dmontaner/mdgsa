##mdGsa.r
##2009-11-16 dmontaner@cipf.es
##2013-03-25 dmontaner@cipf.es


##' @name mdGsa
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases multidimensionalGsa multivariateGsa 
##' 
##' @keywords multidimensional multivariate GSA gene set
##' @seealso \code{\link{uvGsa}}, \code{\link{mdPat}}, \code{glm.fit},
##' \code{p.adjust}
##'
##' @title Multi-Dimensional Gene Set Analysis.
##' 
##' @description
##' Performs a Multi-Variate Gene Set Analysis for two genomic measurements.
##' 
##' @details
##' 'index' must be a numerical matrix or data.frame with at least two columns.
##' 
##' If there are more than three columns,
##' the ranking indexes are taken form the two first one.
##' The remaining columns are used as covariates to correct for
##' within the analysis.
##'
##' Default p-value correction is "BY".
##'
##' In the output data.frame there are three parameters of each type:
##' 'lor', 'pval', \dots
##' one for each of the two genomic conditions analyzed and the third one
##' for the interaction between them.
##'
##' If available, names of the fist two columns of the index matrix
##' are used in the output data.frame.
##' Changing the order of these two first columns will change the report
##' order, but will not change the interpretation of the results.
##' See Montaner et al. (2010) for further details on the algorithm.
##' 
##' @param index ranking index, generally a two column matrix.
##' @param annot an annotation list.
##' @param p.adjust.method p-value adjustment method for multiple testing.
##' @param family see \code{glm}.
##' @param fulltable if TRUE, 'sd', 't' and 'convergence' indicator
##' from the glm fit are included in the output.
##' @param verbose verbose.
##' @param verbosity integer indicating which iterations should be indicated
##' when verbose = TRUE.
##' @param useColnames if TRUE the names of the two first columns of the
##' matrix 'index' are used in the results data.frame.
##' @param \dots further arguments to be pasted to \code{glm.fit},
##' for instance 'weights'.
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
##' Apart from the 'N' coefficient, all other indices appear in triplicate:
##' one coefficient for each genomic condition and a third one for the
##' interaction.
##'
##' @references
##' Montaner et al. (2010)
##' "Multidimensional Gene Set Analysis of Genomic Data."
##' PLoS ONE.
##' 
##' @examples
##' rindexMat <- matrix (rnorm (2000), ncol = 2)
##' colnames (rindexMat) <- c ("genomicVar1", "genomicVar2")
##' rownames (rindexMat) <- paste0 ("gen", 1:1000)
##'
##' annotList <- list (geneSet1 = sample (rownames (rindexMat), size = 10),
##'                    geneSet2 = sample (rownames (rindexMat), size = 15),
##'                    geneSet3 = sample (rownames (rindexMat), size = 20))
##'
##' res <- mdGsa (rindexMat, annotList)
##' res
##' 
##' @export
mdGsa <- function (index, annot, p.adjust.method = "BY",
                   family = quasibinomial(),
                   verbose = TRUE, verbosity = 100,
                   fulltable = FALSE, useColnames = TRUE, ...) {

  
    ## INPUT FORMATTING ########################################################
    
    ##GENE IDS
    if (is.data.frame (index)) {
        index <- as.matrix (index)
    }
    if (is.matrix (index)) {
        genes <- rownames (index)
        rownames (index) <- NULL
    }
    if (is.null (genes)) {
        stop ("no geneIds found. Check rownames in 'index'")
    }
    
    ##BLOCK IDS
    blocks <- names (annot)
    if (is.null (blocks)) {
        stop ("unnamed 'annot'")
    }
    
    
    ## INTERNALS ###############################################################
    
    ##results matrix
    res <- matrix (NA, nrow = length (blocks), ncol = 12+2)
    rownames (res) <- blocks
    colnames (res) <- c("N",
                        "lor.X",   "lor.Y",  "lor.I",
                        "sd.X",     "sd.Y",   "sd.I",
                        "t.X",       "t.Y",    "t.I",
                        "pval.X", "pval.Y", "pval.I",
                        "conv") # "error"
    
    ##store time
    t0 <- proc.time ()
    
    ##set counter
    counter <- 0
    if (verbose) {
        message ("Analyzed blocks:")
    }
    
    ##dimension of index
    index.col.N <- ncol (index)
    if (index.col.N < 2) {
        stop ("'index' must have at least two columns")
    }

  
    ## ALGORITHM ###############################################################

    ## Design matrix:
    ## We use the two first columns of the matrix (as provided by the user)
    ## and include an interaction effect.
    ## Thus, computations are equivalent to glm (Y ~ index[,1] * index[,2])
    X <- cbind (rep (1, times = length (genes)),
                index[,1:2], index[,1] * index[,2])

    ## If extra columns are provided (more than two)
    ## they are understood as co-factors to control for in the model.
    ## Those co-factors are included in the design matrix.
    ## Hence, computations are equivalent to
    ## glm (Y ~ index[,1] * index[,2] + cofactor 1 + cofactor 2 ...
    if (index.col.N > 2) {
        X <- cbind (X, index[,3:index.col.N])  ##covariate correction
    }
    
    colnames (X) <- NULL
    
    for (bl in blocks) {
        B <- as.numeric (genes %in% annot[[bl]])
        res.glm <- glm.fit (x = X, y = B, family = family, ...)
        res.sum <- summary.glm (res.glm)
        res[bl,] <- c (sum (B), res.sum$coefficients[2:4,], res.glm$converged)
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
    ## #########################################################################
  

  
  ##Time
    t1 <- proc.time ()
    if (verbose) {
        message ("time in seconds:")
        print (t1-t0)
    }
    
    ##Convergence
    if (verbose) {
        if (any (res[,"conv"] == 0)) {
            tex <- "The analysis did not converge for some blocks.
                    You may re-run mdGsa using 'fulltable = TRUE' to find them."
            message (gsub ("  +", "", tex))
        }
    }
    
    
  ## p-value ADJUSTMENT #######################################################

    res <- cbind (res,
                 'padj.X' = p.adjust (res[,"pval.X"], method = p.adjust.method),
                 'padj.Y' = p.adjust (res[,"pval.Y"], method = p.adjust.method),
                 'padj.I' = p.adjust (res[,"pval.I"], method = p.adjust.method))
  

    ##OUTPUT ###################################################################
    
    ##reorder columns
    orden <- c ("N",
                "lor.X",   "lor.Y",  "lor.I",
                "pval.X", "pval.Y", "pval.I",
                "padj.X", "padj.Y", "padj.I",
                "sd.X",     "sd.Y",   "sd.I",
                "t.X",       "t.Y",    "t.I",
                "conv")
  
    res <- res[,orden]
  
  ##remove some columns
    if (!fulltable) {
        res <- res[,c ("N",
                       "lor.X",   "lor.Y",  "lor.I",
                       "pval.X", "pval.Y", "pval.I",
                       "padj.X", "padj.Y", "padj.I")]
    }    
    
    ##set the names
    if (useColnames) {
        if (!is.null (colnames (index))) {
            column1 <- colnames (index)[1]
            column2 <- colnames (index)[2]
            colnames (res) <- sub ("X", column1, colnames (res))
            colnames (res) <- sub ("Y", column2, colnames (res))
        }
    }
    
    ##format data.frame
    res <- as.data.frame (res)
    
    ## OUTPUT
    res
}
