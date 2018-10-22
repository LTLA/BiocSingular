# Tests the bpmatrix functionality.
# library(testthat); library(BiocSingular); source("test-bpmatrix.R")

test_that("bpmatrix utilities work correctly", {
    y <- matrix(rnorm(1000), ncol=20)
    colnames(y) <- head(LETTERS, 20)
    x <- BiocSingular:::bpmatrix(y, BiocParallel::SerialParam())
    expect_identical(dim(x), dim(y))
    expect_identical(dimnames(x), dimnames(y))
    expect_identical(length(x), length(y))
})

test_that("bpmatrix multiplication works correctly", {
    y <- matrix(rnorm(1000), ncol=20)
    z <- matrix(rnorm(500), nrow=20)
    ref <- y %*% z

    yp <- BiocSingular:::bpmatrix(y, BiocParallel::SerialParam())
    expect_equal(ref, yp %*% z)

    zp <- BiocSingular:::bpmatrix(z, BiocParallel::SerialParam())
    expect_equal(ref, y %*% zp)

    expect_equal(ref, yp %*% zp)
})

test_that("bpmatrix crossprod works correctly", {
    y <- matrix(rnorm(1000), nrow=20)
    yp <- BiocSingular:::bpmatrix(y, BiocParallel::SerialParam())
    expect_equal(crossprod(y), BiocSingular:::crossprod(yp))

    z <- matrix(rnorm(500), nrow=20)
    ref <- crossprod(y, z)
    expect_equal(ref, BiocSingular:::crossprod(yp, z))

    zp <- BiocSingular:::bpmatrix(z, BiocParallel::SerialParam())
    expect_equal(ref, BiocSingular:::crossprod(y, zp))
    expect_equal(ref, BiocSingular:::crossprod(yp, zp))
})

test_that("bpmatrix tcrossprod works correctly", {
    y <- matrix(rnorm(1000), ncol=20)
    yp <- BiocSingular:::bpmatrix(y, BiocParallel::SerialParam())
    expect_equal(tcrossprod(y), BiocSingular:::tcrossprod(yp))

    z <- matrix(rnorm(500), ncol=20)
    ref <- tcrossprod(y, z)
    expect_equal(ref, BiocSingular:::tcrossprod(yp, z))

    zp <- BiocSingular:::bpmatrix(z, BiocParallel::SerialParam())
    expect_equal(ref, BiocSingular:::tcrossprod(y, zp))
    expect_equal(ref, BiocSingular:::tcrossprod(yp, zp))
})
