##test_getOntology.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_getOntology <- function () {
    require (GO.db)
    system.time (tgt <- Ontology (GOTERM))
    ids <- names (tgt)
    length (ids)
    system.time (cur <- getOntology (ids))
    checkIdentical (tgt, cur)
}

#test_getOntology ()
