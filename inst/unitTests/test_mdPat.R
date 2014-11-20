##test_mdPat.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_mdPat <- function () {
    N <- c (10, 20, 30, 40)
    lor.X <- c (1.45, -0.32, 1.89, -1.66)
    lor.Y <- c (2.36, -1.86, 0.43, -2.01)
    lor.I <- c (0.89, -0.12, 0.24,  3.55)
    pval.X <- c (0.001, 0.002, 0.003, 0.06)
    pval.Y <- c (0.002, 0.003, 0.06,  0.07)
    pval.I <- c (0.003, 0.02,  0.05,  0.08)
    padj.X <- p.adjust (pval.X, "BY")
    padj.Y <- p.adjust (pval.Y, "BY")
    padj.I <- p.adjust (pval.I, "BY")
    
    mdGsa.res <- as.data.frame (cbind (N,
                                       lor.X, lor.Y, lor.I,
                                       pval.X, pval.Y, pval.I,
                                       padj.X, padj.Y, padj.I))
    mdGsa.res

    cur <- mdPat (mdGsa.res)
    
    tgt <- c("q1i", "q3f", "xh", "NS")
    
    checkTrue (all (tgt == cur))
}

#test_mdPat ()
