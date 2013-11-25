##propagateGO.r
##2013-04-03 dmontaner@cipf.es

## To Do
##  Eliminate most of all warnings.

##' @name propagateGO
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords propagate GO gene ontology
## @seealso \code{\link{uvGsa}}, \code{mdPat}, annotMat2list, revList, annotFilter
##'
##' @title Propagate Gene Ontology annotation.
##' 
##' @description
##' Genes annotated under a GO term inherit the annotation for all its ancestors.
##' 
##' @details
##' Uses the library GO.db.
##' 
##' @param annot annotation list or matrix.
##' @param verbose verbose
##'
##' @return A annotation matrix or list with the propagated annotation.
##'
##' @examples
##' mat <- cbind (c("gene1", "gene2"), c("GO:0034390", "GO:0042889"))
##' mat
##' propagateGO (mat)
##'
##' li <- list ('GO:0034390' = "gene1", 'GO:0042889' = "gene2")
##' li
##' propagateGO (li)
##'
##' @importFrom AnnotationDbi as.list
##' @import GO.db
##' 
##' @export
propagateGO <- function (annot, verbose = TRUE) {
  
  if (is.matrix (annot) | is.data.frame (annot)) {
    annot <- propagateGO.matrix (annot, verbose = verbose)
  }
  
  if (is.list (annot)) {
    annot <- annotList2mat (annot)
    annot <- propagateGO.matrix (annot, verbose = verbose)
    annot <- annotMat2list (annot)
  }

  return (annot)
}


################################################################################


## @name propagateGO.matrix
## 
## @author David Montaner \email{dmontaner@@cipf.es}
##
## @keywords propagate GO gene ontology matrix
## 
## @seealso \code{\link{propagateGO}}
##
## @title Propagate Gene Ontology annotations for annotation matrices.
## 
## @description
## Genes annotated under a GO term inherit the annotation for all its ancestors.
## 
## @details
## Uses the library GO.db.
##
## @param annotation annotation list or matrix
## @param verbose verbose.
## 
## @return A annotation matrix with the propagated annotation.
##
## @importFrom AnnotationDbi as.list
## @import GO.db
## @export
propagateGO.matrix <- function (annotation, verbose = TRUE) {

  if (verbose) {
    message ("\n", "Using GO.db version: ", packageDescription ("GO.db", fields = "Version")) #2.3.5
  }
  
  if (verbose) cat ("", fill = TRUE)
  
  t0 <- proc.time ()

  columnas <- colnames (annotation)
  annotation <- as.matrix (annotation)
  dimnames (annotation) <- NULL
  
  duplicados <- duplicated (annotation)
  annotation <- annotation[!duplicados,]
  t1 <- proc.time ()
  if (verbose) cat (c ("   duplicated 1 :", round ((t1 - t0)[1:3], 2)), fill = TRUE)
  
  ancestros <- c (as.list (GOBPANCESTOR), as.list (GOMFANCESTOR), as.list (GOCCANCESTOR))
  t2 <- proc.time ()
  if (verbose) cat (c ("  get ancestors :", round ((t2 - t1)[1:3], 2)), fill = TRUE)
  
  ancestros.ordenados <- ancestros[annotation[,2]]
  t3 <- proc.time ()
  if (verbose) cat (c (" sort ancestors :", round ((t3 - t2)[1:3], 2)), fill = TRUE)
  
  longitudes <- sapply (ancestros.ordenados, length,  USE.NAMES = FALSE)
  t4 <- proc.time ()
  if (verbose) cat (c ("compute lengths :", round ((t4 - t3)[1:3], 2)), fill = TRUE)
  
  heredados <- cbind (rep (annotation[,1], times = longitudes), unlist (ancestros.ordenados))
  rownames (heredados) <- NULL
  t5 <- proc.time ()
  if (verbose) cat (c (" unlist & cbind :", round ((t5 - t4)[1:3], 2)), fill = TRUE)
  
  annotation <- rbind (annotation, heredados)
  t6 <- proc.time ()
  if (verbose) cat (c ("          rbind :", round ((t6 - t5)[1:3], 2)), fill = TRUE)
  
  noees.all <- annotation[,2] != "all"
  annotation <- annotation[noees.all,]
  t7 <- proc.time ()
  if (verbose) cat (c ("     remove all :", round ((t7 - t6)[1:3], 2)), fill = TRUE)
  
  duplicados <- duplicated (annotation)
  annotation <- annotation[!duplicados,]
  t8 <- proc.time ()
  if (verbose) cat (c ("   duplicated 2 :", round ((t8 - t7)[1:3], 2)), fill = TRUE)

  cat ("", fill = TRUE)

  colnames (annotation) <- columnas
  
  return (annotation)
}
