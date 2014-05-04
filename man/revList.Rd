\name{revList}
\alias{revList}
\title{Revert an annotation list.}
\usage{
revList(lis, tag = "listPos")
}
\arguments{
  \item{lis}{annotation list.}

  \item{tag}{substitutes missing list names if any.}
}
\value{
An inverted list.
}
\description{
Inverts a list: names to elements / elements to names
}
\examples{
data (breast)
res <- uvGsa (ranking[,1], annot)
cbind (res, uvPat (res))
}
\author{
David Montaner \email{dmontaner@cipf.es}
}
\seealso{
\code{\link{annotMat2list}}, \code{\link{revList}}
annotMat2list, annotList2mat
}
\keyword{annotation}
\keyword{list}
\keyword{revert}

