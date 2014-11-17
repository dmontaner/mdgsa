##format_annotation.r
##2010-07-16 dmontaner@cipf.es
##2013-03-25 dmontaner@cipf.es


##' @name annotMat2list
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords annotation matrix list
##' @seealso \code{\link{annotList2mat}}, \code{\link{revList}}, \code{split}
##'
##' @title Convert an annotation matrix into an annotation list.
##' 
##' @description
##' Converts an annotation matrix to an annotation list.
##' The annotation matrix should have 2 columns, the first one with the gene ids;
##' the second one with the annotation ids.
##' 
##' @details
##' Each element of the annotation list represents a functional block;
##' it is a character vector containing the gene ids annotated under the function.
##' The names of the list are the annotation ids.
##'
##' @param mat annotation matrix; gene IDs in the first column; block IDs in the second column.
##' 
##' @return An annotation list: elements of the list are vectors of genes;
##' names of the list are Gene Set ids.
##' 
##' @examples
##' mat <- cbind (c("gen1", "gen2", "gen3"), c("Block1", "Block1", "Block2"))
##' annotMat2list (mat)
##' @export
annotMat2list <- function (mat) {
  
  mat <- as.matrix (mat) ## for data.frames
  dimnames (mat) <- NULL
  mat <- unique (mat)
  
  out <- split (x = mat[,1], f = mat[,2], drop = TRUE, sep = "")
  
  ## blocks <- unique (mat[,2])
  ## out <- list ()
  ## for (bk in blocks) {
  ##   out[[bk]] <- mat[mat[,2] == bk, 1]
  ## }
  
  out
}


################################################################################


##' @name annotList2mat
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords annotation matrix list
##' @seealso \code{\link{annotMat2list}}, \code{\link{revList}}, \code{split}
##'
##' @title Convert an annotation list into an annotation matrix.
##'
##' @description
##' Converts an annotation list to an annotation matrix.
##' The annotation matrix should have 2 columns, the first one with the gene ids;
##' the second one with the annotation ids.
##' 
##' @details
##' Each element of the annotation list represents a functional block;
##' it is a character vector containing the gene ids annotated under the function.
##' The names of the list are the annotation ids.
##' 
##' @param lis annotation list.
##' @param tag substitutes missing list names if any.
##' 
##' @return An annotation matrix: the first column contains the gene or feature ids,
##' the second column contains the Gene Set or functional block ids.
##' 
##' @examples
##' lis <- list (Block1 = c("gen1", "gen2"), Block2 = c("gen3"))
##' annotList2mat (lis)
##' 
##' @export
annotList2mat <- function (lis, tag = "listPos") {
  
  ## names (lis): may be NULL or contain NA
  ## tag: for missing list names tagging

  ## list   NAMES  go to column 2
  ## list ELEMENTS go to column 1 
  
  nombres <- names (lis)
  ##
  if (is.null (nombres)) {
    nombres <- paste (tag, 1:length (lis), sep = "")
  } else {
    esna <- is.na (nombres)
    nombres[esna] <- paste (tag, which (esna), sep = "")
  }    

  longitudes <- sapply (lis, length)

  v.nombres <- rep (nombres, times = longitudes)
  v.element <- unlist (lis)

  if (length (v.nombres) != length (v.element)) stop ("Matrix could not be reconstructed. Revise the structure of the input list.")

  salida <- cbind (v.element, v.nombres) ##consistent with annotMat2list
  ##
  dimnames (salida) <- NULL
  salida
}



################################################################################


##' @name revList
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords revert annotation list
##' @seealso \code{\link{annotMat2list}}, \code{\link{revList}} annotMat2list, annotList2mat
##'
##' @title Revert an annotation list.
##'
##' @description
##' Inverts a list: names to elements / elements to names
##' 
##' @param lis annotation list.
##' @param tag substitutes missing list names if any.
##' 
##' @return An inverted list. 
##' 
##' @examples
##' lis <- list (Block1 = c("gen1", "gen2"), Block2 = c("gen1", "gen3"))
##' revList (lis)
##' 
##' @export
revList <- function (lis, tag = "listPos") {
  
  annmat <- annotList2mat (lis)
  annmat <- annmat[,2:1]
  
  salida <- annotMat2list (annmat)
  salida
}
