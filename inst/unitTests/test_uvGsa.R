##test_uvGsa.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_uvGsa <- function () {

    ## rindex is matrix
    rindex <- matrix (rnorm (1000), ncol = 1)
    rownames (rindex) <- paste0 ("gen", 1:1000)
    
    annotList <- list (geneSet1 = sample (rownames (rindex), size = 10),
                       geneSet2 = sample (rownames (rindex), size = 15),
                       geneSet3 = sample (rownames (rindex), size = 20))
    
    res <- uvGsa (rindex, annotList)
    res

    Y <- rownames (rindex) %in% annotList$geneSet1
    lin <- summary (glm (Y ~ rindex, family = quasibinomial ()))

    tgt <- lin$coefficients[-1, c("Estimate", "Pr(>|t|)")]
    cur <- unlist (res[1, 2:3])

    checkEqualsNumeric (tgt, cur)
    

    ## rindex is vector
    rindex <- rnorm (1000)
    names (rindex) <- paste0 ("gen", 1:1000)
    
    annotList <- list (geneSet1 = sample (names (rindex), size = 10),
                       geneSet2 = sample (names (rindex), size = 15),
                       geneSet3 = sample (names (rindex), size = 20))
    
    res <- uvGsa (rindex, annotList)
    res

    Y <- names (rindex) %in% annotList$geneSet1
    lin <- summary (glm (Y ~ rindex, family = quasibinomial ()))

    tgt <- lin$coefficients[-1, c("Estimate", "Pr(>|t|)")]
    cur <- unlist (res[1, 2:3])

    checkEqualsNumeric (tgt, cur)
}

#test_uvGsa ()
