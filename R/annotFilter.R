##annotFilter.r
##2010-07-16 dmontaner@cipf.es
##2013-04-05 dmontaner@cipf.es


##' @name annotFilter
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases
##' 
##' @keywords filter annotation list
##' @seealso \code{\link{annotMat2list}}
##' 
##' @title Checks and filters an annotation list.
##'
##' @description
##' Checks that the annotated genes (those in the annotation list)
##' are consistent with the universe of genes defined by the ranking index.
##' Filters out functional blocks too 'big' or too 'small'.
##' 
##' @details
##' \code{index} is optional.
##' When it is not provided, the annotation lists is just filtered out
##' by the sizes of the blocks of genes defined in the list.
##'
##' If a ranking index is provided its names are assumed to be the universe of genes under study.
##' The genes in the annotation list are compared against those of the ranking index
##' and the ones not belonging to the universe are removed out form the annotation
##' in order to compute the size of each functional block.
##' Then the list is filtered by sizes; too big and too small blocks are removed.
##'
##' No transformation is done over the ranking index or its names (gene IDs).
##'
##' \code{index} may just be a character vector containing the names of the genes in the universe, 
##' that is, the names or rownames of the ranking index.
##' 
##' @param annot an annotation list.
##' @param index ranking index. Vector, matrix or data.frame
##' @param minBlockSize minimum block size kept
##' @param maxBlockSize maximum block size kept
##' @param verbose verbose
##' 
##' @return a filtered annotation list.
##'
##' @examples
##' rindex <- 1:10
##' names (rindex) <- paste ("gene", rindex, sep = "")
##' rindex
##'
##' annot <- list (paste ("gene", 1:2, sep = ""),          ##too small block
##'                paste ("gene", 1:6, sep = ""),          ##right size
##'                paste ("gene", 1:5, sep = ""),          ##too big block
##'                paste ("gene", c(1:3, 1:3), sep = ""))  ##duplicated IDs
##' annot[[2]][1] <- NA
##' annot[[2]][2] <- ""
##' annot[[2]][3] <- "BAD_ID"
##' annot
##'
##' annotFilter (annot, minBlockSize = 3, maxBlockSize = 5)
##' annotFilter (annot, rindex, minBlockSize = 3, maxBlockSize = 5)
##' 
##' @export

annotFilter <- function (annot, index, minBlockSize = 10, maxBlockSize = 500, verbose = TRUE) {
  
  ##BLOCK IDs
  blocks <- names (annot)
  
  if (is.null (blocks))  {
    cat ("\n", "Warning: The annotation list has no names.", sep = "", fill = TRUE)
    cat ("List position will be used instead.", fill = TRUE)
    names (annot) <- paste ("block", 1:length (annot), sep = "_")
    blocks <- names (annot)
  }
  
  pos <- which (is.na (blocks))
  if (length (pos) > 0) {
    cat ("\n", "Warning: Some names are missing for the blocks in the annot list.", sep = "", fill = TRUE)
    cat ("List positions will be used instead.", fill = TRUE)
    names (annot)[pos] <- paste ("block", pos, sep = "_")
    blocks <- names (annot)
  }
  
  if (any (duplicated (blocks))) {
    cat ("\n", "Warning: Some names are duplicated for the blocks in the annot list.", sep = "", fill = TRUE)
  }
  
  ## ###########################################################################

  ##Duplicated genes
  annot.aux <- lapply (annot, unique)
  if (any (sapply (annot, length) != sapply (annot.aux, length))) {
    cat ("\n",  "Warning: Some blocks in the annot list have duplicated genes.", sep = "", fill = TRUE)
    cat ("Duplicated will be removed.", fill = TRUE)
  }
  annot <- annot.aux
  
  ##genes that are NA or ""
  annot.aux <- lapply (annot, setdiff, y = c(NA, ""))
  if (any (sapply (annot, length) != sapply (annot.aux, length))) {
    cat ('\n', 'Warning: Some blocks in the annot gene IDs which are missing or "".', sep = "", fill = TRUE)
    cat ("Those will be removed.", fill = TRUE)
  }
  annot <- annot.aux
  
  ## ###########################################################################
  
  if (!missing (index)) {
    
    if (is.matrix (index) | is.data.frame (index)) {
      gen.universe <- rownames (index)
    } else {
      if (is.character (index)) {
        gen.universe <- index       ##see whether it is convenient this usage of the index variable; may be we should create a gene universe third variable
      } else {
        gen.universe <- names (index) ##assuming that 'index' is a vector
      }
    }
    
    if (is.null (gen.universe)) {
      stop ("Names where not found in index.")
    }
    
    ## #########################################################################
    
    ##duplicated or missing IDs in the Gene Universe
    if (any (is.na (gen.universe)))      cat ('\n', 'Warning: Some genes or feature names form the ranking index have an NA value.',    sep = "", fill = TRUE)
    if (any (gen.universe == ""))        cat ('\n', 'Warning: Some genes or feature names form the ranking index have an "" value.',    sep = "", fill = TRUE)
    if (any (duplicated (gen.universe))) cat ('\n', 'Warning: Some genes or feature names from the ranking index have duplicated IDs.', sep = "", fill = TRUE)
    
    gen.universe <- setdiff (unique (gen.universe), c(NA, ""))
    
    ## #########################################################################
   
    ##Just IDs in the Gene Universe
    genes <- unique (unlist (annot))
    genes.aux <- genes %in% gen.universe
    
    if (any (!genes.aux)) {
      cat ("\n", "Warning: There are genes in the annotation list which are not part of the gene universe defended by the ranking index;", sep = "", fill = TRUE)
      cat ("they will be eliminated form the annotation.", fill = TRUE)
      
      annot.aux <- lapply (annot, intersect, y = gen.universe)  ##Here 'gene.universe' is already unique and has no NA or "" values.
      annot <- annot.aux
    }
    
    ## #########################################################################
    
    ##some stats
    cat ("\n",
         paste (round (100 * sum (gen.universe %in% genes) / length (gen.universe), 2),
                "% of the genes in the index are annotated in the list.", sep = ""),
         sep = "",
         fill = TRUE)
    
  } ##END if missing (index)
  
  
  ##FILTER by size
  annot <- annot.size.filter (annot = annot, minBlockSize = minBlockSize, maxBlockSize = maxBlockSize, verbose = verbose)

  ##OUT
  ##invisible (annot)
  annot
}

################################################################################
################################################################################


## @name annot.size.filter
## @docType 
## @author David Montaner \email{dmontaner@@cipf.es}
## 
## @keywords annotation list filter size
## @seealso \code{\link{annotFilter}}
## 
## @title Filter an annotation list by size.
## 
## @description
## Filter the annotation list by size
##
## @details
## Preferably use annotFilter which makes some more checking.
##
## @param annot an annotation list.
## @param minBlockSize minimum block size
## @param maxBlockSize maximum block size
## @param verbose verbose
##
## @return An annotation list filtered by size.
##
## @export
annot.size.filter <- function (annot, minBlockSize, maxBlockSize, verbose = TRUE) {

  if (minBlockSize < 1) {
    warning ("minBlockSize < 1; Generally is not meaningful to allow for empty blocks")
  }
  
  block.size <- sapply (annot, length)

  is.small <- block.size < minBlockSize
  is.big   <- block.size > maxBlockSize
  
  touse <- !is.big & !is.small 
  annot <- annot[touse]

  if (verbose) {
    cat ("", fill = TRUE)
    cat ("Filtering annotation by size:", fill = TRUE)
    cat (paste ("  ", sum (is.small), "small blocks removed."),            fill = TRUE)
    cat (paste ("  ", sum (is.big),   "big blocks removed."),              fill = TRUE)
    cat (paste ("  ", sum (touse),    "blocks remain in the annotation."), fill = TRUE)
    cat ("", fill = TRUE)
  }

  ###OUT
  annot
}
