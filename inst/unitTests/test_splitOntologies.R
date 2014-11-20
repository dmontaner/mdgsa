##test_splitOntologies.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_splitOntologies <- function () {
    require (GO.db)
    ##help (package = GO.db)

    annot <- list ("GO:0006915" = c ("g1"),
                   "GO:0016020" = c ("g2", "g3"),
                   "GO:0008152" = c ("g1", "g2", "g3"),
                   "GO:0015288" = c ("g4", "g5"))
    annot
    spl <- splitOntologies (annot)
    spl
    
    cur.bp <- names (spl[["bp"]])
    cur.mf <- names (spl[["mf"]])
    cur.cc <- names (spl[["cc"]])
    
    gos <- names (annot)
    onto <- Ontology (GOTERM[gos])
    onto
    
    tgt.bp <- names (onto[onto == "BP"])
    tgt.mf <- names (onto[onto == "MF"])
    tgt.cc <- names (onto[onto == "CC"])

    checkEquals (tgt.bp, cur.bp)
    checkEquals (tgt.mf, cur.mf)
    checkEquals (tgt.cc, cur.cc)
}

#test_splitOntologies ()
