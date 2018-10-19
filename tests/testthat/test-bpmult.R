# Tests bpmult() for parallelized matrix multiplication.
# library(testthat); library(BiocSingular); source("test-bpmult.R")

test_that("bpmult works correctly on two input matrices", {
    # 'y' is bigger than 'x'.
    x <- matrix(runif(1000), ncol=20)
    y <- matrix(runif(5000), nrow=20)
    ref <- x%*%y

    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=MulticoreParam(3)))

    # 'x' is bigger than 'y'.
    x_ <- t(y)
    y_ <- t(x)
    ref <- x_%*%y_

    expect_equal(ref, BiocSingular:::bpmult(x_, y_, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(x_, y_, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult(x_, y_, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[0,]
    ref <- x0 %*% y
    expect_equal(ref, BiocSingular:::bpmult(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult(x0, y, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult(x0, y, BPPARAM=MulticoreParam(3)))

    y0 <- y[,0]
    ref <- x %*% y0
    expect_equal(ref, BiocSingular:::bpmult(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult(x, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult(x, y0, BPPARAM=MulticoreParam(3)))

    x0_ <- x_[0,]
    ref <- x0_%*%y_
    expect_equal(ref, BiocSingular:::bpmult(x0_, y_, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult(x0_, y_, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult(x0_, y_, BPPARAM=MulticoreParam(3)))

    y0_ <- y_[,0]
    ref <- x_%*%y0_
    expect_equal(ref, BiocSingular:::bpmult(x_, y0_, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult(x_, y0_, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult(x_, y0_, BPPARAM=MulticoreParam(3)))

    # Both are empty.
    x0 <- x[,0]
    y0 <- y[0,]
    ref <- x0 %*% y0
    expect_equal(ref, BiocSingular:::bpmult(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult(x0, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult(x0, y0, BPPARAM=MulticoreParam(3)))
})

test_that("bpmult works with vectors", {
    x <- matrix(runif(1000), ncol=20)
    y <- runif(20)
    ref <- x%*%y

    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=MulticoreParam(3)))

    # Flipping it.
    x_ <- y
    y_ <- t(x)
    ref <- x_%*%y_

    expect_equal(ref, BiocSingular:::bpmult(x_, y_, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(x_, y_, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult(x_, y_, BPPARAM=MulticoreParam(3)))

    # Two vectors.
    z <- runif(20)
    ref <- z %*% y
    expect_equal(ref, BiocSingular:::bpmult(z, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(z, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult(z, y, BPPARAM=MulticoreParam(3)))
})

