##test_uvPat.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_uvPat <- function () {
    uvGsa.res <- as.data.frame (list (N    = c (10, 20, 30, 40),
                                      lor  = c (1.45, -0.32, 1.89, -1.66),
                                      pval = c (0.001, 0.002, 0.05, 0.06)))
    uvGsa.res[,"padj"] <- p.adjust (uvGsa.res$pval, "BY")
    uvGsa.res[,"pat"] <- uvPat (uvGsa.res)
    uvGsa.res

    cur <- uvGsa.res[,"pat"]    

    tgt <- sign (uvGsa.res[,"lor"]) * as.numeric (uvGsa.res[,"padj"] < 0.05)

    checkEqualsNumeric (tgt, cur)
}

#test_uvPat ()
