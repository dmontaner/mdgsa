##plotMdGsa.r
##2009-11-16 dmontaner@cipf.es
##2013-09-25 dmontaner@cipf.es


##' @name plotMdGsa
## @docType 
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
## @aliases 
##' 
##' @keywords multidimensional multivariate GSA plot
##' @seealso \code{\link{mdGsa}}, \code{\link{mdPat}}, \code{\link{mdPat}}, \code{ellipsoidPoints}
##'
##' @title Plot Multi-Dimensional Gene Set
##' 
##' @description
##' Plots confidence region for a Gene Set in a two dimensional space.
##' 
##' @details
##' Black dots show all genes in the dataset. Red stars show the genes in the Gene Set.
##' Blue axis show the "center" of the distribution of all genes;
##' blue ellipse shows the confidence region for all genes.
##' Red axis show the "center" of the distribution of the genes in the Gene Set;
##' red ellipse shows the confidence region for genes in the Gene Set.
##' 
##' @param index matrix or data frame with the two columns of ranking statistics.
##' @param block matrix or data frame with \strong{gene} ids in the first column and \strong{gene set} ids in the second column.
##' @param cr level of the confidence region.
##' @param pch plotting character for all genes.
##' @param pch.block plotting character for the genes in the gene set or functional block.
##' @param lwd line width. Used when drawing ellipses and other lines. 
##' @param col.all color used to represent all genes.
##' @param col.block color used to represent the genes in the gene set being plotted.
##' @param project if TRUE projection over the axis are displayed for the genes of the gene set.
##' @param col.proj color used to plot the projection.
##' @param diagonals if TRUE diagonals are plotted.
##' @param col.diag color used to plot the diagonals.
##' @param \dots arguments to be passed to plot
##' 
##' @return A plot.
##' 
## @examples
##' 
## @import cluster
##' @importFrom cluster ellipsoidPoints
##' 
##' 
##' @export

plotMdGsa <- function (index, block, cr = 0.95,
                       pch = ".", pch.block = 20, lwd = 2, 
                       col.all = "blue", col.block = "red",
                       project = FALSE, col.proj = "green",
                       diagonals = FALSE, col.diag = "grey",
                       ...) {
  
  genes <- rownames (index)
  
  ##all genes
  plot (index[,1], index[,2], xlab = colnames (index)[1], ylab = colnames (index)[2], pch = pch, lwd = lwd, ...)
  
  ##genes in the Gene Set
  points (index[block, 1], index[block, 2], col = col.block, pch = pch.block, lwd = lwd, ...)
  
  ##ellipses; needs the library (cluster)
  C.ls <- cov      (index[block,])
  m.ls <- colMeans (index[block,])
  d2 <- qchisq(cr, df = 2)
  lines (ellipsoidPoints(C.ls, d2, loc = m.ls), col = col.block, lwd = lwd, ...)
  abline (v = m.ls[1], h = m.ls[2], col = col.block, lwd = lwd, ...)
  ##
  C.ls <- cov (index)
  m.ls <- colMeans (index)
  d2 <- qchisq(cr, df = 2)
  lines (ellipsoidPoints(C.ls, d2, loc = m.ls), col = col.all, lwd = lwd, ...)
  abline (v = m.ls[1], h = m.ls[2], col = col.all, lwd = lwd, ...)
  
  ##projections
  if (project) {
    points (index[block, 1], rep (m.ls[2], times = length (block)), col = col.proj, pch = "|", lwd = lwd, ...)
    points (rep (m.ls[1], times = length (block)), index[block, 2], col = col.proj, pch = "_", lwd = lwd, ...)
  }
  
  ##diagonals
  if (diagonals) {
    abline (a = m.ls[2]-m.ls[1], b =  1, col = col.diag, lwd = lwd, ...)
    abline (a = m.ls[2]+m.ls[1], b = -1, col = col.diag, lwd = lwd, ...)
  }
}

################################################################################
