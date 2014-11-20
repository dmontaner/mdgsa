##test_get_GO_KEGG_names.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_getGOnames <- function () {
    require (GO.db)
    system.time (tgt <- Term (GOTERM))
    ids <- names (tgt)
    system.time (cur <- getGOnames (ids))
    checkIdentical (tgt, cur)
}

test_getKEGGnames <- function () {
    require (KEGG.db)
    system.time (tgt <- unlist (as.list (KEGGPATHID2NAME)))
    ids <- names (tgt)
    system.time (cur <- getKEGGnames (ids))
    checkIdentical (tgt, cur)
}

#test_getGOnames ()
#test_getKEGGnames ()
