# Tests runIrlba().
# library(testthat); library(BiocSingular); source("setup.R"); source("test-irlba.R")

library(irlba)
library(BiocParallel)
old <- bpparam()
register(SerialParam())

set.seed(9000)
test_that("IRLBA works on input matrices", {

    y <- matrix(rnorm(50000), ncol=200)
    set.seed(100)
    out <- runIrlba(y, k=5, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5)
    expect_equal_svd(out, ref[c("d", "u", "v")])

    # Handles truncation.
    set.seed(100)
    out <- runIrlba(y, k=10, nv=5, nu=3, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=10, nu=3)
    ref$v <- ref$v[,1:5]
    expect_equal_svd(out, ref[c("d", "u", "v")])

    set.seed(100)
    out <- runIrlba(y, k=1, nv=5, nu=3, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5, nu=3)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref[c("d", "u", "v")])

    # Handles zeroes.
    out0 <- runIrlba(y, k=0, nv=0, nu=0)
    expect_equal(out0$d, numeric(0))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(out0$v, ref$v[,0])
})

set.seed(9001)
test_that("IRLBA works on thin matrices", {
    y <- matrix(rnorm(10000), ncol=10)
    set.seed(200)
    out <- runIrlba(y, k=3, fold=1)
    set.seed(200)
    ref <- irlba(y, nv=3, nu=3)
    expect_equal_svd(out, ref)

    # Handles truncation.
    set.seed(200)
    out <- runIrlba(y, k=5, nv=3, nu=2, fold=1)
    set.seed(200)
    ref <- irlba(y, nu=2, nv=5)
    ref$v <- ref$v[,1:3]
    expect_equal_svd(out, ref)

    set.seed(200)
    out <- runIrlba(y, k=1, nv=3, nu=2, fold=1)
    set.seed(200)
    ref <- irlba(y, nu=2, nv=3)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)
})

set.seed(9002)
test_that("IRLBA works on fat matrices", {
    y <- matrix(rnorm(10000), nrow=10)
    set.seed(300)
    out <- runIrlba(y, k=4, fold=1)
    set.seed(300)
    ref <- irlba(y, nu=4, nv=4)
    expect_equal_svd(out, ref)

    # Handles truncation.
    set.seed(300)
    out <- runIrlba(y, k=5, nv=4, nu=2, fold=1)
    set.seed(300)
    ref <- irlba(y, nu=2, nv=5)
    ref$v <- ref$v[,1:4]
    expect_equal_svd(out, ref)

    set.seed(300)
    out <- runIrlba(y, k=1, nv=6, nu=2, fold=1)
    set.seed(300)
    ref <- irlba(y, nu=2, nv=6)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)
})

set.seed(9003)
test_that("IRLBA works with alternative multiplication", {
    # Lower tol to avoid numeric precision issues.
    y <- matrix(rnorm(50000), ncol=100)
    set.seed(100)
    out <- runIrlba(y, k=5, BPPARAM=MulticoreParam(2), tol=1e-8, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5, tol=1e-8)
    expect_equal_svd(out, ref[c("d", "u", "v")], tol=1e-6)

    # Handles truncation.
    set.seed(100)
    out <- runIrlba(y, k=5, nv=2, nu=3, BPPARAM=MulticoreParam(2), tol=1e-8, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5, nu=3, tol=1e-8)
    ref$v <- ref$v[,1:2]
    expect_equal_svd(out, ref[c("d", "u", "v")], tol=1e-6)
})

set.seed(9004)
test_that("IRLBA works with centering and scaling", {
    y <- matrix(rnorm(10000), ncol=50)
    center <- runif(ncol(y))
    scale <- NULL # runif(ncol(y)) # irlba seems to be broken with scaling.

    set.seed(100)
    ref <- irlba(y, k=5, center=center, scale=scale, tol=1e-8)
    set.seed(100)
    out <- runIrlba(y, k=5, center=center, scale=scale, fold=Inf, tol=1e-8) 
    expect_equal_svd(out, ref)

    ry <- scale(y, center=center, scale=FALSE)
    set.seed(100)
    ref2 <- irlba(ry, nv=5, nu=5, tol=1e-8)
    expect_equal_svd(out, ref2)

    # Works with the cross-product.
    y <- matrix(rnorm(10000), ncol=20)
    center <- runif(ncol(y))
    scale <- NULL #runif(ncol(y))
    
    ry <- scale(y, center=center, scale=FALSE)
    set.seed(200)
    ref <- irlba(ry, nu=5, nv=5, tol=1e-8)
    set.seed(200)
    out <- runIrlba(y, k=5, center=center, scale=scale, fold=1, tol=1e-8)
    expect_equal_svd(out, ref, tol=1e-6)

    # Works with the alternative multiplication.
    set.seed(100)
    out <- runIrlba(ry, k=5, BPPARAM=MulticoreParam(2), tol=1e-8)
    set.seed(100)
    ref <- irlba(y, nv=5, center=center, scale=scale, tol=1e-8)
    expect_equal_svd(out, ref, tol=1e-6)
})

register(old)
