##test_mdGsa.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_mdGsa <- function () {

    rindexMat <- matrix (rnorm (2000), ncol = 2)
    colnames (rindexMat) <- c ("genomicVar1", "genomicVar2")
    rownames (rindexMat) <- paste0 ("gen", 1:1000)
    
    annotList <- list (geneSet1 = sample (rownames (rindexMat), size = 10),
                       geneSet2 = sample (rownames (rindexMat), size = 15),
                       geneSet3 = sample (rownames (rindexMat), size = 20))
    
    res <- mdGsa (rindexMat, annotList)
    res

    Y <- rownames (rindexMat) %in% annotList$geneSet1
    lin <- summary (glm (Y ~ rindexMat[,1] * rindexMat[,2], family = quasibinomial ()))

    tgt <- lin$coefficients[-1, c("Estimate", "Pr(>|t|)")]
    cur <- unlist (res[1, 2:7])

    checkEqualsNumeric (tgt, cur)
}

#test_mdGsa ()
