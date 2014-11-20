##test_pval2index.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)

test_pval2index <- function () {
    my.statistic <- rnorm (1000)
    my.pvalue <- 2 * pnorm (my.statistic)
    my.pvalue[my.pvalue > 1] <- 2 - my.pvalue[my.pvalue > 1]
    
    cur <- pval2index (pval = my.pvalue, sign = my.statistic)

    tgt <- (-1) * log (my.pvalue) * sign (my.statistic)

    checkEqualsNumeric (tgt, cur)
}

#test_pval2index ()
