##test_annot_format.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_annotMat2list <- function () {
    mat <- cbind (c("gen1", "gen2", "gen3"), c("Block1", "Block1", "Block2"))
    lis <- annotMat2list (mat)
    ##
    checkTrue (is.list (lis))
    checkIdentical (length (lis), length (unique (mat[,2])))
}

test_annotList2mat <- function () {
    lis <- list (Block1 = c("gen1", "gen2"), Block2 = c("gen3"))
    mat <- annotList2mat (lis)
    ##
    checkTrue (is.matrix (mat))
    checkIdentical (length (lis), length (unique (mat[,2])))
}

test_revList <- function () {
    lis <- list (Block1 = c("gen1", "gen2"), Block2 = c("gen3"))
    mat <- annotList2mat (lis)
    rev <- revList (lis)
    ##
    checkTrue (is.matrix (mat))
    checkTrue (is.list (rev))
    checkIdentical (length (lis), length (unique (mat[,2])))
    checkIdentical (length (rev), length (unique (mat[,1])))
}

## test_annotMat2list ()
## test_annotList2mat ()
## test_revList ()

