% Generated by roxygen2 (4.0.2): do not edit by hand
\name{annotMat2list}
\alias{annotMat2list}
\title{Convert an annotation matrix into an annotation list.}
\usage{
annotMat2list(mat)
}
\arguments{
\item{mat}{annotation matrix; gene IDs in the first column;
block IDs in the second column.}
}
\value{
An annotation list: elements of the list are vectors of genes;
names of the list are Gene Set ids.
}
\description{
Converts an annotation matrix to an annotation list.
The annotation matrix should have 2 columns, the first one with the
gene ids; the second one with the annotation ids.
}
\details{
Each element of the annotation list represents a functional block;
it is a character vector containing the gene ids annotated under the
function.
The names of the list are the annotation ids.
}
\examples{
mat <- cbind (c("gen1", "gen2", "gen3"), c("Block1", "Block1", "Block2"))
annotMat2list (mat)
}
\author{
David Montaner \email{dmontaner@cipf.es}
}
\seealso{
\code{\link{annotList2mat}}, \code{\link{revList}}, \code{split}
}
\keyword{annotation}
\keyword{list}
\keyword{matrix}

