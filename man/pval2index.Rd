\name{pval2index}
\alias{pval2index}
\title{Transform p-values in into a ranking index.}
\usage{
  pval2index(pval, sign, log = TRUE, verbose = TRUE)
}
\arguments{
  \item{pval}{a vector or matrix of p-values.}

  \item{sign}{a vector or matrix of signs associated to the
  p-values.}

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
  By default the transformation is (-1) * log (pval) * sign
  (sign). When log = FALSE the transformation is (1 - pval)
  * sign (sign).
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

