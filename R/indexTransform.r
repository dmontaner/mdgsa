##indexTransform.r
##2013-03-25 dmontaner@cipf.es

## To Do:
##  index.tdist ???

index.standardize <- function (index) {
  m <- mean (index, na.rm = TRUE)
  s <- sd   (index, na.rm = TRUE)
  index <- (index - m) / s
  return (index)
}

index.normalize <- function (index) {
  ##Standardize to normal distribution (ie. quantiles of a normal distribution)

  if (!is.numeric (index)) stop ("index is not numeric")  #check that index is numeric
    
  res <- qqnorm (index, plot.it = FALSE)$x
  if (!is.null (names (index))){
    names (res) <- names (index)
  }
  return (res)
}

################################################################################

##' @name indexTransform
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords transform ranking index
## @seealso \code{\link{uvGsa}}, \code{mdPat}, annotMat2list, revList, annotFilter
##'
##' @title Transform ranking index distribution.
##' 
##' @description
##' The function performs a transformation of the ranking index so that its distribution
##' is suitable as independent variable of the logistic regression model.
##' 
##' @details
##' Works for vector, matrices and data.frames.
##' In the case of matrices the function transforms each column separately from the other ones.
##'
##' Two methods are currently implemented:
##' \describe{
##'   \item{- normalize:}{transforms the index into quantiles of a normal distribution.}
##'   \item{- standardize:}{performs an statistical standardization by subtracting
##'   the mean and dividing by the standard deviation.}
##' }
##'
##' @param index a ranking index; a numeric vector, matrix or data.frame.
##' @param method transformation method. See details.
##' 
##' @return A transformed index. Its class will be that of the the input object.
##'
##' @examples
##' myIndex <- runif (1000)
##' myTransformedIndex <- indexTransform (myIndex)
##' plot (myIndex, myTransformedIndex)
##' 
##' @export
##' 
indexTransform <- function (index, method = "normalize") {
  
  metodos <- c("normalize", "standardize", "none")
  if (!method %in% metodos) {
    stop (paste ("method should be one of:", paste (metodos, collapse = ", ")))
  }
  
  if (method == "normalize") {
    if (is.matrix (index) | is.data.frame (index)) {
      index <- apply (index, 2, index.normalize)
    } else {
      index <- index.normalize (index)
    }
  }
  
  if (method == "standardize") {
    if (is.matrix (index) | is.data.frame (index)) {
      index <- apply (index, 2, index.standardize)
    } else {
      index <- index.standardize (index)
    }
  }
  
  ## if (method == "none") {
  ##   no change
  ## }
  
  ##return
  return (index)
}
