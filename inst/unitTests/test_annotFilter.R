##test_annotFilter.r
##2014-11-19 dmontaner@cipf.es
##UnitTesting

## rm (list = ls ())
## R.version.string ##"R version 3.1.1 (2014-07-10)"
## library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.1"
## library (RUnit); packageDescription ("RUnit", fields = "Version") #"0.4.27"
## ## help (package = RUnit)


test_annotFilter <- function () {
    rindex <- 1:10
    names (rindex) <- paste ("gene", rindex, sep = "")
    ##
    annot <- list (paste ("gene", 1:2, sep = ""),          ##too small block
                   paste ("gene", 1:6, sep = ""),          ##right size
                   paste ("gene", 1:5, sep = ""),          ##too big block
                   paste ("gene", c(1:3, 1:3), sep = ""))  ##duplicated IDs
    annot[[2]][1] <- NA
    annot[[2]][2] <- ""
    annot[[2]][3] <- "BAD_ID"
    ##
    an1 <- annotFilter (annot, minBlockSize = 3, maxBlockSize = 5)
    an2 <- annotFilter (annot, rindex, minBlockSize = 3, maxBlockSize = 5)
    
    lon1 <- sapply (an1, length)
    lon2 <- sapply (an2, length)
    
    checkTrue (all (2 < lon1 & lon1 < 6))
    checkTrue (all (2 < lon2 & lon2 < 6))
    
    checkTrue (all (unlist (an2) %in% names (rindex)))
}

#test_annotFilter ()
