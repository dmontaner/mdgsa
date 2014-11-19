##get_GO_KEGG_names.r
##2013-04-03 dmontaner@cipf.es
##2013-09-26 dmontaner@cipf.es


##' @name getGOnames
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords GO ontology names
##' @seealso \code{\link{propagateGO}}, \code{\link{goLeaves}},
##' \code{\link{splitOntologies}}, \code{\link{getKEGGnames}},
##' \code{\link{getOntology}}
##' 
##' @title Get Gene Ontology names
##'
##' @description
##' Finds the GO name form GO id.
##' 
##' @details
##' Uses the library GO.db.
##'
##' \code{x} may be a \code{data.frame}.
##' In such case, GO ids are expected in its row names.
##' 
##' @param x a character vector of GO ids.
##' @param verbose verbose.
##'
##' @return A character vector with the corresponding GO names.
##'
##' @examples
##' getGOnames (c("GO:0000018", "GO:0000038", "BAD_GO"))
##' 
##' @import DBI
##' @import GO.db
##'
##' @export

getGOnames <- function (x, verbose = TRUE) {
    
    if (is.data.frame (x) | is.matrix (x)) {
        if (verbose) message ("Using row names of the input matrix.")
        x <- rownames (x)
    }
    
    if (verbose) {
        message ("Using GO.db version: ",
                 packageDescription ("GO.db", fields = "Version")) #2.9.0
    }
    
    ##go id to ontology
    micon <- GO_dbconn ()
    tabla <- dbReadTable (micon, "go_term")
    tabla <- tabla[,c("go_id", "term")]
    tabla <- tabla[tabla$go_id != "all",]
    
    id2name <- tabla[,"term"]
    names (id2name) <- tabla[,"go_id"]
    
    ##my go ids
    res <- id2name[x]
    
    if (any (is.na (res))) {
        warning (sum (is.na (res)),
                 " GOids where not found; missing names generated.")
    }
    
    res
}


################################################################################
################################################################################


##' @name getKEGGnames
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords KEGG names
##' @seealso \code{\link{getGOnames}}
##' 
##' @title Get KEGG  names
##' 
##' @description
##' Finds the KEGG name form KEGG id.
##' 
##' @details
##' Uses the library KEGG.db.
##'
##' \code{x} may be a \code{data.frame}.
##' In such case, GO ids are expected in its row names.
##' 
##' @param x a character vector of KEGG ids.
##' @param verbose verbose.
##'
##' @return A character vector with the corresponding KEGG names.
##'
##' @examples
##' getKEGGnames (c("00010", "00020", "BAD_KEGG"))
##' 
##' @import DBI
##' @import KEGG.db
##'
##' @export

getKEGGnames <- function (x, verbose = TRUE) {

    if (is.data.frame (x) | is.matrix (x)) {
        if (verbose) message ("Using row names of the input matrix.")
        x <- rownames (x)
    }
    
    if (verbose) {
        message ("Using KEGG.db version: ",
                 packageDescription ("KEGG.db", fields = "Version")) #2.9.0
    }

    ##kegg id to kegg name
    micon <- KEGG_dbconn ()
    tabla <- dbReadTable (micon, "pathway2name")
    ##tabla <- tabla[,c("path_id", "path_name")]
    ##anything to filter out?
    
    id2name <- tabla[,"path_name"]
    names (id2name) <- tabla[,"path_id"]
    
    ##my kegg ids
    res <- id2name[x]
    
    if (any (is.na (res))) {
        warning (sum (is.na (res)),
                 " KEEGids where not found; missing names generated.")
    }
    
    res
}
