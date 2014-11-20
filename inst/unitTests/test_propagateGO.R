##test_propagateGO.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_propagateGO <- function () {
    require (GO.db)
    mat <- cbind (c("gene1", "gene2"), c("GO:0034390", "GO:0042889"))
    cur <- propagateGO (mat)[,2]
    cur <- unique (cur)
    
    tgt <- c(c("GO:0034390", "GO:0042889"),
             setdiff (unlist (as.list (GOBPANCESTOR[c("GO:0034390", "GO:0042889")])), "all"))
    tgt <- unique (tgt)
    
    checkEquals (tgt, cur)
}

#test_propagateGO ()
