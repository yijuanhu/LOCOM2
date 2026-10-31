context("Testing `locom2` function")
library(LOCOM2)
library(testthat)

data("throat.otu.table.filter")
data("throat.meta.filter")
data("throat.otu.taxonomy")


Y <- ifelse(throat.meta.filter$SmokingStatus == "NonSmoker", 0, 1)
C <- ifelse(throat.meta.filter$Sex == "Male", 0, 1)

# test
test_that("`locom2` function provides expected results", {
    system.time(res <- locom2(formula = throat.otu.table.filter | C ~ Y, seed = 123, n.perm.max=100))
    expect_equivalent(length(res$detected.otu.Wald), 6)
})
