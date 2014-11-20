##test_indexTransform.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## help (package = RUnit)

test_indexTransform <- function () {
    rindex <- rnorm (100)
    tgt <- order (rindex)
    cur <- order (indexTransform (rindex))
    checkIdentical (tgt, cur)
}

#test_indexTransform ()
