# LOCOM2

A logistic regression model for testing differential abundance in compositional microbiome data

The LOCOM2 package includes the `locom2' function that allows you to test
(1). whether any OTU (or taxon) is associated with the trait of interest with FDR control, based on log ratios of relative abundances between pairs of taxa, and
(2). whether the whole community is associated with the trait (a global test), based on the harmonic mean method for combining individual p-values
The tests accommodate continuous, discrete (binary or categorical), and multivariate traits, and allow adjustment for confounders.

LOCOM2 extends LOCOM (Hu et al., 2022, PNAS) in the following ways:
-. accommodating both relative abundance and read count data for OTUs;
-. refining the weighting scheme in LOCOM to eliminate confounding by library size;
-. incorporating a series of adjustments to ensure stable and reliable inference, even under extreme conditions such as rare taxa and highly unbalanced case–control designs;
-. replacing the computationally intensive permutation procedure with a Wald-type test (using a fixed 1,000 permutation replicates). 
For more details see He et al. (2026) <doi:10.64898/2026.04.07.716976>.

This package also includes the 'Com2seq' function that combines two sequencing datasets (e.g., from 16S marker-gene and shotgun metagenomic sequencing) for integrative analysis of differential abundance. For more details see Yue et al. (2025) <doi.org/10.1080/01621459.2025.2516205>.


## Installation

```r
devtools::install_github("yijuanhu/LOCOM2")
```

## An example of using the locom2 function

Apply 'locom2' to a dataset from the study of smoking effects on the microbiome of the human upper respiratory tract (Charlson et al., 2010):

```r
library(LOCOM2)
data("throat.otu.table.filter")
data("throat.meta.filter")
data("throat.otu.taxonomy")

Y <- ifelse(throat.meta.filter$SmokingStatus == "NonSmoker", 0, 1)
C <- ifelse(throat.meta.filter$Sex == "Male", 0, 1)
meta <- cbind(Y, C)

# running locom2
res <- locom2(formula = throat.otu.table.filter | C ~ Y, data = meta, seed = 1, n.cores = 4)
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

## An example of using the com2seq function

Apply 'com2seq' to a simulated dataset:

```r
library(LOCOM2)
data("otu.table.1")
data("otu.table.2")
data("meta.1")
data("meta.2")


# restricting to a subset of samples and OTUs for speed
sub.sam.1 <- c(31:60)
sub.sam.2 <- c(1:20, 51:60)
otu.table.sub.1 <- otu.table.1[sub.sam.1, 1:70]
otu.table.sub.2 <- otu.table.2[sub.sam.2, 1:70]
meta.sub.1 <- meta.1[sub.sam.1,]
meta.sub.2 <- meta.2[sub.sam.2,]

# running com2seq
res.com2seq <- com2seq(formula.1 = otu.table.sub.1 | C ~ as.factor(Y), data.1 = meta.sub.1, 
                       formula.2 = otu.table.sub.2 | C ~ as.factor(Y), data.2 = meta.sub.2,
                       seed = 123, n.cores = 4, n.perm.max = 1000)
res.com2seq$p.global.omni 
res.com2seq$detected.taxa.omni 
```