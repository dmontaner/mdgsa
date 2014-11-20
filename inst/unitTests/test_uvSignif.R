##test_uvSignif.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_uvSignif <- function () {
    uvGsa.res <- as.data.frame (list (N    = c (10, 20, 30, 40),
                                      lor  = c (1.45, -0.32, 1.89, -1.66),
                                      pval = c (0.001, 0.002, 0.05, 0.06)))
    uvGsa.res[,"padj"] <- p.adjust (c (0.001, 0.002, 0.05, 0.06), "BY")
    ##uvGsa.res[,"pat"] <- uvPat (uvGsa.res)
    uvGsa.res

    uvPat (uvGsa.res)
    
    tgt <- uvGsa.res[1:2,]
    cur <- uvSignif (uvGsa.res)

    checkEquals (tgt, cur)
}

#test_uvSignif ()
