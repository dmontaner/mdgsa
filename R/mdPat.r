##mdPat.r
##2009-11-16 dmontaner@cipf.es
##2013-03-27 dmontaner@cipf.es


##' @name mdPat
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords multidimensional multivariate GSA pattern
##' @seealso \code{\link{mdGsa}}, \code{\link{uvPat}}
##'
##' @title Multi-Dimensional Gene Set Analysis Pattern Classification.
##' 
##' @description
##' Classifies significant patterns form a Multi-Variate Gene Set Analysis.
##' 
##' @details
##' Sign of the three 'lor' and p-values are used to classify functional blocks.
##' The classification is done in the two dimensional space previously analyzed
##' by mdGsa.
##' 
##' @param gsaout data.frame; output from mdGsa.
##' @param cutoff p-value cutoff for considering significant a Gene Set.
##' @param pvalue p-value column to be used. Default is named "padj" as in mdGsa output.
##' 
##' @return A character vector indicating the pattern associated to each Gene Set.
##' 
##' @references \href{fhttp://www.plosone.org/article/info\%3Adoi\%2F10.1371\%2Fjournal.pone.0010348}{Montaner et al. 2010 "Multidimensional Gene Set Analysis of Genomic Data." PLoS ONE.}
##' 
##' @examples
##' data (breast)
##' res <- mdGsa (ranking, annot)
##' cbind (res, mdPat (res))
##' 
##' @export
mdPat <- function (gsaout, cutoff = 0.05, pvalue = "padj") {
  
  ##mdGsa Patterns
  patron <- c ( 1,  1,  1, "q1i",
                1,  0,  1, "q1i",
                0,  1,  1, "q1i",
               -1, -1,  1, "q3i",
               -1,  0,  1, "q3i",
                0, -1,  1, "q3i",
               -1,  1, -1, "q2i",
               -1,  0, -1, "q2i",
                0,  1, -1, "q2i",
                1, -1, -1, "q4i",
                1,  0, -1, "q4i",
                0, -1, -1, "q4i",
                0,  0,  1, "b13",
                0,  0, -1, "b24", #corrected from the paper
                1,  1,  0, "q1f",
               -1, -1,  0, "q3f",
               -1,  1,  0, "q2f",
                1, -1,  0, "q4f",
                1,  0,  0, "xh",
               -1,  0,  0, "xl",
                0,  1,  0, "yh",
                0, -1,  0, "yl",
                0,  0,  0, "NS",
                1,  1, -1, "q1f", #not in the paper
               -1, -1, -1, "q3f", #not in the paper
               -1,  1,  1, "q2f", #not in the paper
                1, -1,  1, "q4f") #not in the paper
  
  patron <- matrix (patron, ncol = 4, byrow = TRUE)
  ##colnames (patron) <- c("alpha", "beta", "gamma", "pattern")
  mdGsaPattern <- patron[,4]
  mdGsaPattern.name <- apply (patron[,1:3], 1, paste, collapse = "") 
  names (mdGsaPattern) <- mdGsaPattern.name

  
  ##COLUMNS
  lors <- sub ("lor", "", grep ("lor.", colnames (gsaout), value = TRUE))
  pvals <- sub (pvalue, "", grep (paste (pvalue, ".", sep = ""), colnames (gsaout), value = TRUE))

  tags <- intersect (lors, pvals)

  if (length (tags) != 3) {
    stop ("lor and pvalue columns could not be matched")
  } else {
    lors <- paste ("lor", tags, sep = "")
    pvals <- paste (pvalue, tags, sep = "")
  }

  ## Log Odds Ratio sign * {1, 0} for significant and no significant
  significance <- sign (gsaout[,lors]) * as.numeric (gsaout[,pvals] < cutoff)
  colnames (significance) <- c("alpha", "beta", "gamma")

  
  ##pattern assignment
  significance.id <- apply (significance, 1, paste, collapse = "") 
  res <- mdGsaPattern[significance.id]
  names (res) <- rownames (significance)

  ##non relevant patterns
  res[is.na (res)] <- "NR"  #CHECK THIS
  
  ##return
  return (res)
}
