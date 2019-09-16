# Tests dispatch in the runSVD() function.
# library(BiocSingular); library(testthat); source("test-svd.R")

set.seed(99001)
test_that("exact method dispatches correctly", {
    y <- matrix(rnorm(10000), ncol=100, nrow=100)
    expect_identical(runExactSVD(y, k=3, fold=Inf), runSVD(y, k=3, BSPARAM=ExactParam(fold=Inf)))
    expect_identical(runExactSVD(y, k=4, fold=1), runSVD(y, k=4, BSPARAM=ExactParam(fold=1)))
    expect_identical(runExactSVD(y, k=5, deferred=TRUE), runSVD(y, k=5, BSPARAM=ExactParam(deferred=TRUE)))
    
    # Other arguments are passed along.    
    expect_identical(runExactSVD(y, k=5, center=TRUE), runSVD(y, k=5, center=TRUE, BSPARAM=ExactParam()))
    expect_identical(runExactSVD(y, k=5, scale=TRUE), runSVD(y, k=5, scale=TRUE, BSPARAM=ExactParam()))

    # Missing method dispatches to exact.
    expect_identical(runExactSVD(y, k=6), runSVD(y, k=6))
})

set.seed(99002)
test_that("irlba method dispatches correctly", {
    y <- matrix(rnorm(10000), ncol=100, nrow=100)
    
    set.seed(100); left <- runIrlbaSVD(y, k=3, fold=Inf)
    set.seed(100); right<- runSVD(y, k=3, BSPARAM=IrlbaParam(fold=Inf))
    expect_identical(left, right)

    set.seed(200); left <- runIrlbaSVD(y, k=4, fold=1)
    set.seed(200); right <- runSVD(y, k=4, BSPARAM=IrlbaParam(fold=1))
    expect_identical(left, right)

    set.seed(300); left <- runIrlbaSVD(y, k=5, deferred=TRUE)
    set.seed(300); right <- runSVD(y, k=5, BSPARAM=IrlbaParam(deferred=TRUE))
    expect_identical(left, right)
    
    # Other arguments are passed along.    
    set.seed(400); left <- runIrlbaSVD(y, k=5, center=TRUE)
    set.seed(400); right <- runSVD(y, k=5, center=TRUE, BSPARAM=IrlbaParam())
    expect_identical(left, right)

# irlba doesn't handle scale=TRUE, for reasons unknown to me.
#    set.seed(500); left <- runIrlbaSVD(y, k=5, scale=TRUE)
#    set.seed(500); right <- runSVD(y, k=5, scale=TRUE, BSPARAM=IrlbaParam())
#    expect_identical(left, right)
})

set.seed(99003)
test_that("random method dispatches correctly", {
    y <- matrix(rnorm(10000), ncol=100, nrow=100)
    
    set.seed(100); left <- runRandomSVD(y, k=3, fold=Inf)
    set.seed(100); right<- runSVD(y, k=3, BSPARAM=RandomParam(fold=Inf))
    expect_identical(left, right)

    set.seed(200); left <- runRandomSVD(y, k=4, fold=1)
    set.seed(200); right <- runSVD(y, k=4, BSPARAM=RandomParam(fold=1))
    expect_identical(left, right)

    set.seed(300); left <- runRandomSVD(y, k=5, deferred=TRUE)
    set.seed(300); right <- runSVD(y, k=5, BSPARAM=RandomParam(deferred=TRUE))
    expect_identical(left, right)
    
    # Other arguments are passed along.    
    set.seed(400); left <- runRandomSVD(y, k=5, center=TRUE)
    set.seed(400); right <- runSVD(y, k=5, center=TRUE, BSPARAM=RandomParam())
    expect_identical(left, right)

    set.seed(500); left <- runRandomSVD(y, k=5, scale=TRUE)
    set.seed(500); right <- runSVD(y, k=5, scale=TRUE, BSPARAM=RandomParam())
    expect_identical(left, right)
})

set.seed(99003)
test_that("fast automatic method dispatches correctly", {
    y <- matrix(rnorm(10000), ncol=100, nrow=100)

    set.seed(100); left <- runSVD(y, k=3, BSPARAM=FastAutoParam())
    set.seed(100); right<- runSVD(y, k=3, BSPARAM=IrlbaParam())
    expect_identical(left, right)
    
    set.seed(200); left <- runSVD(y, k=3, BSPARAM=FastAutoParam(fold=1))
    set.seed(200); right<- runSVD(y, k=3, BSPARAM=IrlbaParam(fold=1))
    expect_identical(left, right)

    dy <- DelayedArray(y)
    set.seed(300); left <- runSVD(dy, k=3, BSPARAM=FastAutoParam())
    set.seed(300); right<- runSVD(dy, k=3, BSPARAM=RandomParam())
    expect_identical(left, right)

    dfm <- DeferredMatrix(y)
    set.seed(400); left <- runSVD(dfm, k=3, BSPARAM=FastAutoParam())
    set.seed(400); right<- runSVD(dfm, k=3, BSPARAM=IrlbaParam())
    expect_identical(left, right)
})
