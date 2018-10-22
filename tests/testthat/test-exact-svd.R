# Tests runExactSVD().
# library(testthat); library(BiocSingular); source("setup.R"); source("test-exact-svd.R")

test_that("exact SVD works on square matrices", {
    y <- matrix(rnorm(10000), ncol=100, nrow=100)
    out <- runExactSVD(y, fold=Inf)
    ref <- svd(y)
    expect_equal(out, ref)

    # Handles truncation.
    out <- runExactSVD(y, k=10, nv=5, nu=3, fold=Inf)
    ref <- svd(y, nv=5, nu=3)
    ref$d <- ref$d[1:10]
    expect_equal(out, ref)

    # Handles zeroes.
    out0 <- runExactSVD(y, k=0, nv=0, nu=0)
    expect_equal(out0$d, numeric(0))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(out0$v, ref$v[,0])
})

test_that("exact SVD works on thin matrices", {
    y <- matrix(rnorm(10000), ncol=10)
    out <- runExactSVD(y, fold=1)
    ref <- svd(y)
    expect_equal_svd(out, ref)

    # Handles truncation.
    out <- runExactSVD(y, k=5, nv=6, nu=2, fold=1)
    ref <- svd(y, nu=2, nv=6)
    ref$d <- ref$d[1:5]
    expect_equal_svd(out, ref)

    # Handles zeroes.
    out0 <- runExactSVD(y, k=0, nv=0, nu=0, fold=1)
    expect_equal(out0$d, numeric(0))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(out0$v, ref$v[,0])
})

test_that("exact SVD works on fat matrices", {
    y <- matrix(rnorm(10000), nrow=10)
    out <- runExactSVD(y, fold=1)
    ref <- svd(y)
    expect_equal_svd(out, ref)

    # Handles truncation.
    out <- runExactSVD(y, k=5, nv=6, nu=2, fold=1)
    ref <- svd(y, nu=2, nv=6)
    ref$d <- ref$d[1:5]
    expect_equal_svd(out, ref)

    # Handles zeroes.
    out0 <- runExactSVD(y, k=0, nv=0, nu=0, fold=1)
    expect_equal(out0$d, numeric(0))
    expect_equivalent(out0$u, ref$u[,0])
    expect_equivalent(out0$v, ref$v[,0])
})

test_that("exact SVD works with centering and scaling", {
    y <- matrix(rnorm(10000), ncol=50)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))

    ry <- scale(y, center=center, scale=scale)
    ref <- svd(ry)
    out <- runExactSVD(y, center=center, scale=scale, fold=Inf)
    expect_equal_svd(out, ref)

    # Works with the crossproduct.
    y <- matrix(rnorm(10000), ncol=10)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))

    ry <- scale(y, center=center, scale=scale)
    ref <- svd(ry)
    out <- runExactSVD(y, center=center, scale=scale, fold=1)
    expect_equal_svd(out, ref)
})
