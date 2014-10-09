##goLeaves.r
##2013-04-04 dmontaner@cipf.es.com
##2013-09-26 dmontaner@cipf.es.com

##Leaf nodes. One or more of the nodes in the GO
##DAG will have no children. These are called leaves or leaf nodes.

##' @name goLeaves
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords leaves child go terms
##' @seealso \code{\link{uvGsa}}, \code{\link{uvPat}}, \code{\link{propagateGO}}, \code{\link{pval2index}}
##'
##' @title Keep just leaf nodes from the Gene Ontology DAG.
##'
##' @description
##' Cuts significant terms and filters out all redundant GO terms
##' from a list of \code{uvGsa} results.
##' 
##' @details
##' Uses the library GO.db to find the 'ancestors' of each GO term.
##' Those ancestors are discarded form the \code{uvGsa} results.
##'
##' Alternativeliy, the function may also take a character vector of GO ids
##' in the \code{gsaout} parameter.
##' In such case the function returns also a character vector of GO ids,
##' containing just the GO terms being "leaves" of the original set.
##' 
##' @param gsaout data.frame; output from uvGsa.
##' @param cutoff p-value cutoff for considering significant a Gene Set.
##' @param pvalue p-value column to be used. Default is named "padj" as in uvGsa output.
##' @param statistic name of the column containing the log odds ratio from the uvGsa analysis.
##' @param sort if TRUE the output data.frame is ordered according to significance.
##' @param verbose verbose
##' 
##' @return The input data.frame but keeping just the 'significant' and 'non redundant' GO terms.
##' 
##' @importFrom AnnotationDbi as.list
##' @import GO.db
##'
##' @export

goLeaves <- function (gsaout, cutoff = 0.05, pvalue = "padj", statistic = "lor", verbose = TRUE, sort = TRUE) {
    
    if (verbose) {
        message ("\n", "Using GO.db version: ", packageDescription ("GO.db", fields = "Version")) #2.3.5
    }
    
    ancestros <- c (as.list (GOBPANCESTOR), as.list (GOMFANCESTOR), as.list (GOCCANCESTOR))  
  
    ## ###########################################################################

    ## real gsaout
    if (is.data.frame (gsaout) | is.matrix (gsaout)) {
  
        pat <- uvPat (gsaout, cutoff = cutoff, pvalue = pvalue, statistic = statistic)
  
        touse <- pat %in% c (-1, 1)
        gsaout <- gsaout[touse, , drop = FALSE]
        pat <- pat[touse]
  
        is.up <- pat == 1
        is.do <- pat ==-1
  
        gos.up <- rownames (gsaout)[is.up]
        gos.do <- rownames (gsaout)[is.do]
  
        anc.up <- unique (unlist (ancestros[gos.up]))
        anc.do <- unique (unlist (ancestros[gos.do]))

        leaves.up <- setdiff (gos.up, anc.up)
        leaves.do <- setdiff (gos.do, anc.do)

        leaves <- c(leaves.up, leaves.do)
        touse <- rownames (gsaout) %in% leaves
        gsaout <- gsaout[touse, , drop = FALSE]

        ##sort
        if (sort) {
            myindex <- pval2index (pval = gsaout[,pvalue], sign = gsaout[,statistic], log = TRUE, verbose = verbose)
            orden <- order (myindex, decreasing = TRUE)
            gsaout <- gsaout[orden,]
        }
    }

    ## ###########################################################################

    ## character verctor of GO ids
    if (is.character (gsaout)) {
        anc <- unique (unlist (ancestros[gsaout]))
        leaves <- setdiff (gsaout, anc)
        gsaout <- leaves
    }
  
    ## RETURN
    return (gsaout)
}
