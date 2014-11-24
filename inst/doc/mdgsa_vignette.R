## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----, include = FALSE---------------------------------------------------
library (knitr) ## this is needed when R CMD build
opts_chunk$set (message = FALSE)

## ------------------------------------------------------------------------
library (ALL)
data (ALL)

des.mat <- model.matrix (~ 0 + mol.biol, data = ALL)
colnames (des.mat) <- c("ALL", "BCR", "E2A", "NEG", "NUP", "p15")
head (des.mat)

## ------------------------------------------------------------------------
library (limma)
cont.mat <- makeContrasts (ALL-NEG, BCR-NEG, levels = des.mat)
cont.mat

fit <- lmFit (ALL, design = des.mat)
fit <- contrasts.fit (fit, cont.mat)
fit <- eBayes (fit)

## ------------------------------------------------------------------------
fit$t[1:3,]
fit$p.value[1:3,]

## ------------------------------------------------------------------------
library (hgu95av2.db)
anmat <- toTable (hgu95av2PATH)
anmat[1:3,]

## ------------------------------------------------------------------------
fit$t[1:3, "BCR - NEG"]
fit$p.value[1:3, "BCR - NEG"]

## ------------------------------------------------------------------------
library (mdgsa)

## ------------------------------------------------------------------------
rindex <- pval2index (pval = fit$p.value[,"BCR - NEG"], sign = fit$t[,"BCR - NEG"])
rindex <- indexTransform (rindex)
rindex[1:3]

## ------------------------------------------------------------------------
plot (fit$t[,"BCR - NEG"], rindex)

## ------------------------------------------------------------------------
plot (fit$p.value[,"BCR - NEG"], rindex)

## ------------------------------------------------------------------------
anmat[1:3,]
annot <- annotMat2list (anmat)
length (annot)

## ------------------------------------------------------------------------
lapply (annot[1:3], head, n= 3)

## ------------------------------------------------------------------------
annot <- annotFilter (annot, rindex)

## ----uvGsa---------------------------------------------------------------
res.uv <- uvGsa (rindex, annot)

## ------------------------------------------------------------------------
res.uv[1:3,]

## ------------------------------------------------------------------------
res.uv[,"pat"] <- uvPat (res.uv, cutoff = 0.05)
table (res.uv[,"pat"])

## ----, results = "hide"--------------------------------------------------
res.uv[,"KEGG"] <- getKEGGnames (res.uv)

## ------------------------------------------------------------------------
res <- uvSignif (res.uv)
res[,c("pat", "KEGG")]

## ------------------------------------------------------------------------
fit$t[1:3,]
fit$p.value[1:3,]

## ------------------------------------------------------------------------
rindex <- pval2index (pval = fit$p.value, sign = fit$t)
rindex <- indexTransform (rindex)
rindex[1:3,]

## ------------------------------------------------------------------------
plot (rindex, pch = ".")

## ----mdGsa---------------------------------------------------------------
res.md <- mdGsa (rindex, annot)

## ------------------------------------------------------------------------
res.md[1:3,]

## ------------------------------------------------------------------------
res.md[,"pat"] <- mdPat (res.md)
table (res.md[,"pat"])

## ----, results = "hide"--------------------------------------------------
res.md[,"KEGG"] <- getKEGGnames (res.md)
res.md[1:3,]

## ------------------------------------------------------------------------
Q3 <- rownames (res.md)[res.md$pat == "q3f"]
Q3
plotMdGsa (rindex, block = annot[[Q3]], main = res.md[Q3, "KEGG"])

## ------------------------------------------------------------------------
BI <- rownames (res.md)[res.md$pat == "b13"]
plotMdGsa (rindex, block = annot[[BI]], main = res.md[BI, "KEGG"])

## ------------------------------------------------------------------------
rownames (res.md)[res.md$pat == "yl"]
YL <- "00280"
plotMdGsa (rindex, block = annot[[YL]], main = res.md[YL, "KEGG"])

## ------------------------------------------------------------------------
Q2 <- rownames (res.md)[res.md$pat == "q2f"]
Q2
plotMdGsa (rindex, block = annot[[Q2]], main = res.md[Q2, "KEGG"])

## ------------------------------------------------------------------------
rownames (res.md)[res.md$pat == "yh"]
YH <- "05332"
plotMdGsa (rindex, block = annot[[YH]], main = res.md[YH, "KEGG"])

## ------------------------------------------------------------------------
sessionInfo()

## ----, echo = FALSE, results = "asis"------------------------------------
toLatex (sessionInfo())

