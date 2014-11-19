##transferIndex.r
##2013-10-10 dmontaner@cipf.es


##' @name transferIndex
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords transfer index to genes miRNAs
##' @seealso \code{\link{uvGsa}}, \code{\link{mdGsa}},
##' \code{\link{indexTransform}}
##'
##' @title Transfer a ranking index from regulatory elements,
##' such as miRNAs, to genes.
##' 
##' @description
##' Transfers the ranking index information from some regulatory elements
##' to their target genes.
##' It can be for instance used to infer gene regulation levels from miRNA
##' differential expression levels.
##' Afterwords, the inferred gene index can be explored in a
##' univariate gene set analysis using \code{uvGsa},
##' or in a multivariate gene set analysis using \code{mdGsa}.
##' 
##' @details
##' The function was in principle designed to transfer a ranking index
##' defined for miRNAs to their target genes,
##' but it may also be used to deal with some other regulatory elements
##' as transcription factors for instance.
##'
##' If we have for instance a t statistic accounting for miRNAs
##' differential expression levels,
##' we may add up (or average) all t values from the miRNAs that
##' regulate a gene in order to derive an interference score for that gene.
##' Thus we may get an interference index
##' which may be interpreted in terms of functional blocks or gene sets.
##'
##' \code{index} may be a matrix.
##' In such case just the first column will be used.
##'
##' \code{targets} should be a named list.
##' The names of the list have to coincide (or at least overlap)
##' with the names in the \code{index}.
##' The elements of the list are character vectors containing gene IDs.
##' In the miRNA example, each element of \code{targets}
##' will contain the genes regulated by a miRNA,
##' and the miRNA IDs will appear in the names of the \code{targets}
##' list as well as in the ranking index.
##'
##' @param index ranking index, generally a numerical named vector.
##' @param targets a list describing the regulation blocks. See details.
##' @param method the method to collapse the gene information.
##' Currently 'sum' and 'average' are available.
##' @param transferMatrix if true the transference matrix is returned
##' instead of the ranking statistic.
##' @param verbose verbose.
##' 
##' @return A ranking index transferred to the IDs in the elements of the
##' list \code{targets}.
##'
##' @import Matrix
##'
##' @examples
##'
##' ## miRNA to gene list (targets)
##' targets <- list (mirna1 = "g1", mirna2 = c("g1", "g2"), mirna3 = "g3")
##' 
##' ## original index
##' index <- rnorm (5)
##' names (index) <- paste0 ("mirna", 1:5)
##' 
##' ##transfered index
##' tindex <- transferIndex (index, targets)
##' 
##' ##transformed (normalized) index
##' rindex <- indexTransform (tindex)
##' 
##' ##NOTICE in this case:
##' index["mirna1"] + index["mirna2"] == tindex["g1"]
##' 
##' \dontrun{
##' res <- uvgsa (rindex, annot)
##' }
##' 
##' @export

transferIndex <- function (index, targets, method = "sum",
                           verbose = TRUE, transferMatrix = FALSE) {
  
    ## CHECK
    if (!method %in% c("sum", "average")) {
        stop (paste ('method should be: "sum" or "averaged"'))
    }
    
    ## INPUT FORMATTING
    if (is.list (targets)) {
        genes <- sort (unique (unlist (targets)))
        mirnas.in.targets <- names (targets)
    } else {
        stop ("'targets' has to be a list.")
    }  
    
    if (is.data.frame (index)) {
        index <- as.matrix (index)
    }
    if (is.matrix (index)) {
        mirnas.in.index <- rownames (index)
        index <- index[,1]                  ## use just the first column
    }
    if (is.vector (index)) {
        mirnas.in.index <- names (index)
    }
    if (is.null (mirnas.in.index)) {
        stop ("no IDs found. Check names or rownames in 'index'")
    }
    
    ## miRNA universe
    mirnas <- intersect (mirnas.in.index, mirnas.in.targets)
    mirnas.just.index <- setdiff (mirnas.in.index, mirnas.in.targets)
    mirnas.just.targets <- setdiff (mirnas.in.targets, mirnas.in.index)
    
    if (verbose) {
        message ("  ", length (mirnas), "miRNAs with annotated targets")
        message ("  ", length (mirnas.just.index), "miRNAs without targets")
        message ("  ", length (mirnas.just.targets),
                 "miRNAs with targets but not in the ranking index")
    }
    
    
    ## #########################################################################
    ## TRANSFER: size and speed should be fine for 20.000 genes and 5.000 miRNAs
    ## #########################################################################
    
    ## Matrix
    mat <- Matrix (0, nrow = length (genes), ncol = length (mirnas))
    rownames (mat) <- genes
    colnames (mat) <- mirnas
    ##
    for (mi in mirnas) {
        mat[,mi] <- genes %in% targets[[mi]] 
    }
    
    ## transfer the index to each gene
    mat <- t(t (mat) * index[mirnas])
    
    ## collapse genes
    if (method == "sum") {
        tindex <- rowSums (mat)
        names (tindex) <- rownames (mat)
        ## rowSums method for Matrix does not keep names as
        ## the default method does.
    }
    if (method == "average") {
        tindex <- rowMeans (mat)
        names (tindex) <- rownames (mat)
    }
    
    
    ## OUTPUT
    if (transferMatrix) {
        mat
    } else {
        tindex
    }
}
