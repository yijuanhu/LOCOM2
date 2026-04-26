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
    res <- locom2(otu.table = throat.otu.table.filter, Y = Y, C = C, seed = 123)
    expect_equivalent(length(res$detected.otu.Wald), 5)
})
