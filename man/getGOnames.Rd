\name{getGOnames}
\alias{getGOnames}
\title{Get Gene Ontology names}
\usage{
  getGOnames(x, verbose = TRUE)
}
\arguments{
  \item{x}{a character vector of GO ids.}

  \item{verbose}{verbose.}
}
\value{
  A character vector with the corresponding GO names.
}
\description{
  Finds the GO name form GO id.
}
\details{
  Uses the library GO.db.
}
\examples{
getGOnames (c("GO:0000018", "GO:0000038", "BAD_GO"))
}
\author{
  David Montaner \email{dmontaner@cipf.es}
}
\seealso{
  \code{\link{propagateGO}}, \code{\link{goLeaves}},
  \code{\link{splitOntologies}}, \code{\link{getKEGGnames}}
}
\keyword{GO}
\keyword{names}
\keyword{ontology}

