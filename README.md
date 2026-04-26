# LOCOM2

A logistic regression model for testing differential abundance in compositional microbiome data

The LOCOM2 package allows you to test
(1). whether any OTU (or taxon) is associated with the trait of interest with FDR control, based on log ratios of relative abundances between pairs of taxa, and
(2). whether the whole community is associated with the trait (a global test), based on the harmonic mean method for combining individual p-values
The tests accommodate continuous, discrete (binary or categorical), and multivariate traits, and allow adjustment for confounders.

LOCOM2 extends LOCOM (Hu et al., 2022, PNAS) in the following ways:
-. accommodating both relative abundance and read count data for OTUs;
-. refining the weighting scheme in LOCOM to eliminate confounding by library size;
-. incorporating a series of adjustments to ensure stable and reliable inference, even under extreme conditions such as rare taxa and highly unbalanced case–control designs;
-. replacing the computationally intensive permutation procedure with a Wald-type test (using a fixed 1,000 permutation replicates).


## Installation

```r
devtools::install_github("yijuanhu/LOCOM2")
```

## An example of using the locom function

Apply 'locom2' to a dataset from the study of smoking effects on the microbiome of the human upper respiratory tract (Charlson et al., 2010):

```r
library(LOCOM2)
data("throat.otu.table.filter")
data("throat.meta.filter")
data("throat.otu.taxonomy")

Y <- ifelse(throat.meta.filter$SmokingStatus == "NonSmoker", 0, 1)
C <- ifelse(throat.meta.filter$Sex == "Male", 0, 1)

# running locom2
res <- locom2(otu.table = throat.otu.table.filter, Y = Y, C = C, seed = 1)
length(res$p.otu.Wald)
length(res$detected.otu.Wald) 
res$detected.otu.Wald 
res$p.otu.Wald[res$detected.otu.Wald] 
res$p.global

# Summarizing results
w <- match(res$detected.otu.Wald, names(res$p.otu.Wald))
o <- w[order(res$p.otu.Wald[w])] # o <- order(res$p.otu.Wald)[1:10] # top 10 otus

summary.tab.locom2 <- data.frame(p.value = signif(res$p.otu.Wald[o], 3),
                                 q.value = signif(res$q.otu.Wald[o], 3),
                                 mean.freq = signif(colMeans(throat.otu.table.filter/rowSums(throat.otu.table.filter))[o], 3),
                                 prop.presence = signif(colMeans(throat.otu.table.filter > 0)[o], 3),
                                 beta = signif(res$beta[o], 3), 
                                 beta.error = signif(res$beta.var[o], 3),
                                 otu.num = names(res$p.otu.Wald)[o],
                                 otu.name = throat.otu.taxonomy[as.numeric(names(res$p.otu.Wald)[o]) + 1],
                                 row.names = NULL)
summary.tab.locom2
```
