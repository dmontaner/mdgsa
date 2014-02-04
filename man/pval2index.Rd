\name{pval2index}
\alias{pval2index}
\title{Transform p-values in into a ranking index.}
\usage{
  pval2index(pval, sign, names = NULL, log = TRUE,
    verbose = TRUE)
}
\arguments{
  \item{pval}{a vector or matrix of p-values.}

  \item{sign}{a vector or matrix of signs associated to the
  p-values.}

  \item{names}{a character vector of the names of the
  features}

  \item{log}{= TRUE}

  \item{verbose}{verbose}
}
\value{
  A transformed index. A vector or matrix, depending on the
  input parameters.
}
\description{
  After a genomic test, p-values are numerical indexes
  which account for certain biological characteristic. By
  definition p-values are bounded between zero and one, but
  this may not be suitable as an index. Moreover, p-values
  are always derived form a statistic which sign may be
  important. The function helps transforming the p-value
  and its associated statistic into a ranking index.
}
\details{
  The default transformation is (-1) * log (pval) * sign
  (sign). When log = FALSE the transformation is (1 - pval)
  * sign (sign).

  Missing values are allowed and return NA values. Some
  infinite values may be returned if pval has values
  exactly equal to zero.

  By default the names of the output vector (or row names
  in a matrix) are those of \code{pval} or \code{sign}. If
  \code{names} is provided, then it is used instead.
}
\examples{
my.statistic <- rnorm (1000)
my.pvalue <- 2 * pnorm (my.statistic)
my.pvalue[my.pvalue > 1] <- 2 - my.pvalue[my.pvalue > 1]

index <- pval2index (pval = my.pvalue, sign = my.statistic)

#par (mfrow = c (1,2))
#plot (my.statistic, my.pvalue)
#plot (my.statistic, index)
}
\author{
  David Montaner \email{dmontaner@cipf.es}
}
\seealso{
  \code{\link{indexTransform}}
}
\keyword{index}
\keyword{p-value}
\keyword{ranking}

