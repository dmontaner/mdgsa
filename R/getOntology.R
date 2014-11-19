##getOntology.r
##2014-07-07 dmontaner@cipf.es


##' @name getOntology
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords GO ontology names
##' @seealso \code{\link{getGOnames}}, \code{\link{propagateGO}},
##' \code{\link{goLeaves}}, \code{\link{splitOntologies}},
##' \code{\link{getKEGGnames}}
##' 
##' @title Get GO term Ontology 
##'
##' @description
##' Finds the ontology of a term from its id.
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
##' getOntology (c("GO:0000018", "GO:0005788", "BAD_GO"))
##' 
##' @import DBI
##' @import GO.db
##'
##' @export

getOntology <- function (x, verbose = TRUE) {
    
    if (is.data.frame (x) | is.matrix (x)) {
        if (verbose) message ("Using rownames of x")
        x <- rownames (x)
    }
    
    if (verbose) {
        message ("\n", "Using GO.db version: ",
                 packageDescription ("GO.db", fields = "Version")) #2.9.0
    }
    
    ##go id to ontology
    micon <- GO_dbconn ()
    tabla <- dbReadTable (micon, "go_term")
    tabla <- tabla[,c("go_id", "ontology")]
    tabla <- tabla[tabla$ontology != "universal",]
    
    if2ontology <- tabla[,"ontology"]
    names (if2ontology) <- tabla[,"go_id"]
    
    ##my go ids
    res <- if2ontology[x]
    
    if (any (is.na (res))) {
        warning (sum (is.na (res)),
                 " GOids where not found; missing ontologies generated.")
    }

    ## OUTPUT
    res
}
