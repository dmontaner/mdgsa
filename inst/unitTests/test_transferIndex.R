##test_transferIndex.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)
## ## help (package = mdgsa)

test_transferIndex <- function () {
    targets <- list (mirna1 = "g1", mirna2 = c ("g1", "g2"), mirna3 = "g3")
    
    index <- rnorm (5)
    names (index) <- paste0 ("mirna", 1:5)

    cur <- transferIndex (index, targets)
    
    rev <- revList (targets)

    tgt1 <- sum (index[rev[["g1"]]])
    
    checkEqualsNumeric (tgt1, cur[1])
}

#test_transferIndex ()
