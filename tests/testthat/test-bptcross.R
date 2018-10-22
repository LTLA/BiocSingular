# Tests bptcross() for parallelized transposed cross products.
# library(testthat); library(BiocSingular); source("test-bptcross.R")

test_that("bpcross_x_by_row works correctly", {
    x <- matrix(runif(5000), ncol=20)
    y <- matrix(runif(1000), ncol=20)
    ref <- tcrossprod(x, y)

    expect_equal(ref, BiocSingular:::bptcross_x_by_row(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross_x_by_row(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bptcross_x_by_row(x, y, BPPARAM=MulticoreParam(3)))

    # Only one matrix.
    ref <- tcrossprod(x)
    expect_equal(ref, BiocSingular:::bptcross_x_by_row(x, y=NULL, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross_x_by_row(x, y=NULL, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bptcross_x_by_row(x, y=NULL, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[,0]
    ref <- tcrossprod(x0)
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x0, y=NULL, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x0, y=NULL, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    ref <- tcrossprod(x0, y)
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x0, y, BPPARAM=MulticoreParam(2)))

    y0 <- y[0,]
    ref <- tcrossprod(x, y0)
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x, y0, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    y0 <- y[0,]
    ref <- tcrossprod(x0, y0)
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_x_by_row(x0, y0, BPPARAM=MulticoreParam(2)))
})

test_that("bpcross_y_by_col works correctly", {
    x <- matrix(runif(5000), ncol=20)
    y <- matrix(runif(1000), ncol=20)
    ref <- tcrossprod(x, y)

    expect_equal(ref, BiocSingular:::bptcross_y_by_row(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross_y_by_row(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bptcross_y_by_row(x, y, BPPARAM=MulticoreParam(3)))

    # Handles empties okay.
    x0 <- x[0,]
    ref <- tcrossprod(x0, y)
    expect_equivalent(ref, BiocSingular:::bptcross_y_by_row(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_y_by_row(x0, y, BPPARAM=MulticoreParam(2)))

    y0 <- y[0,]
    ref <- tcrossprod(x, y0)
    expect_equivalent(ref, BiocSingular:::bptcross_y_by_row(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_y_by_row(x, y0, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    y0 <- y[0,]
    ref <- tcrossprod(x0, y0)
    expect_equivalent(ref, BiocSingular:::bptcross_y_by_row(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_y_by_row(x0, y0, BPPARAM=MulticoreParam(2)))
})

test_that("bptcross_by_col works correctly", {
    x <- matrix(runif(1000), ncol=20)
    y <- matrix(runif(5000), ncol=20)
    ref <- tcrossprod(x, y)

    expect_equal(ref, BiocSingular:::bptcross_by_col(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross_by_col(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bptcross_by_col(x, y, BPPARAM=MulticoreParam(3)))

    # Only one matrix.
    ref <- tcrossprod(x)
    expect_equal(ref, BiocSingular:::bptcross_by_col(x, y=NULL, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross_by_col(x, y=NULL, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bptcross_by_col(x, y=NULL, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[,0]
    ref <- tcrossprod(x0)
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x0, y=NULL, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x0, y=NULL, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    ref <- tcrossprod(x0, y)
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x0, y, BPPARAM=MulticoreParam(2)))

    y0 <- y[0,]
    ref <- tcrossprod(x, y0)
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x, y0, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    y0 <- y[0,]
    ref <- tcrossprod(x0, y0)
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bptcross_by_col(x0, y0, BPPARAM=MulticoreParam(2)))
})

test_that("bptcross overall function works correctly", {
    x <- matrix(runif(300), ncol=15)
    ref <- tcrossprod(x)
    expect_equal(ref, BiocSingular:::bptcross(x, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross(x, BPPARAM=MulticoreParam(2))) 
    expect_equal(ref, BiocSingular:::bptcross(x, BPPARAM=MulticoreParam(3)))

    x <- matrix(runif(150), ncol=15)
    y <- matrix(runif(300), ncol=15)
    ref <- tcrossprod(x, y)
    expect_equal(ref, BiocSingular:::bptcross(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bptcross(x, y, BPPARAM=MulticoreParam(3)))

    # Handles vector inputs.
    xv <- runif(20)
    ref <- tcrossprod(xv)
    expect_equal(ref, BiocSingular:::bptcross(xv, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross(xv, BPPARAM=MulticoreParam(2))) 

    y <- matrix(runif(300), ncol=20)
    ref <- tcrossprod(xv, y)
    expect_equal(ref, BiocSingular:::bptcross(xv, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bptcross(xv, y, BPPARAM=MulticoreParam(2)))

    x <- matrix(runif(150), ncol=15)
    yv <- runif(10)
    expect_error(ref <- tcrossprod(x, yv), "non-conformable arguments")
    expect_error(BiocSingular:::bptcross(x, yv, BPPARAM=MulticoreParam(2)), "non-conformable arguments")
})
