% Generated by roxygen2 (4.0.2): do not edit by hand
\name{uvGsa}
\alias{univariateGsa}
\alias{uvGsa}
\title{Uni-Variate Gene Set Analysis.}
\usage{
uvGsa(index, annot, p.adjust.method = "BY", family = quasibinomial(),
  verbose = TRUE, verbosity = 100, fulltable = FALSE, ...)
}
\arguments{
\item{index}{ranking index, generally a numerical named vector.}

\item{annot}{an annotation list.}

\item{p.adjust.method}{p-value adjustment method for multiple testing.}

\item{family}{see \code{glm.fit}.}

\item{verbose}{verbose.}

\item{verbosity}{integer indicating which iterations should be indicated
if verbose = TRUE.}

\item{fulltable}{if TRUE, 'sd', 't' and 'convergence'
indicator from the glm fit are included in the output.}

\item{\dots}{further arguments to be pasted to \code{glm.fit},
for instance 'weights'.}
}
\value{
A \code{data.frame} with a row for each Gene Set or block.
Columns are:
\describe{
  \item{\code{N}:}{number of genes annotated to the Gene Set.}
  \item{\code{lor}:}{log Odds Ratio estimated for the Gene Set.}
  \item{\code{pval}:}{p-values associated to each log Odds Ratio.}
  \item{\code{padj}:}{adjusted p-values.}
  \item{\code{sd}:}{standard deviations associated to each log Odds Ratio.}
  \item{\code{t}:}{t statistic associated to each log Odds Ratio.}
}
}
\description{
Performs a Uni-Variate Gene Set Analysis using a logistic regression model.
}
\details{
'index' may also be a numerical \code{matrix} or \code{data.frame}.
If such a matrix has more than one column,
the ranking index is taken form the first one.
The remaining columns are used as covariates to correct for
within the analysis.

Default p-value correction is "BY".
}
\examples{
rindex <- rnorm (1000)
names (rindex) <- paste0 ("gen", 1:1000)

annotList <- list (geneSet1 = sample (names (rindex), size = 10),
                   geneSet2 = sample (names (rindex), size = 15),
                   geneSet3 = sample (names (rindex), size = 20))

res <- uvGsa (rindex, annotList)
res
}
\author{
David Montaner \email{dmontaner@cipf.es}
}
\seealso{
\code{\link{mdGsa}}, \code{\link{uvPat}},
\code{glm.fit}, \code{p.adjust}
}
\keyword{GSA}
\keyword{gene}
\keyword{set}
\keyword{univariate}

