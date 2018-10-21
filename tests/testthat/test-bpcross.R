# Tests bpcross() for parallelized cross products.
# library(testthat); library(BiocSingular); source("test-bpcross.R")

test_that("bpcross_x_by_col works correctly", {
    x <- matrix(runif(1000), nrow=20)
    y <- matrix(runif(5000), nrow=20)
    ref <- crossprod(x, y)

    expect_equal(ref, BiocSingular:::bpcross_x_by_col(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross_x_by_col(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpcross_x_by_col(x, y, BPPARAM=MulticoreParam(3)))

    # Only one matrix.
    ref <- crossprod(x)
    expect_equal(ref, BiocSingular:::bpcross_x_by_col(x, y=NULL, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross_x_by_col(x, y=NULL, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpcross_x_by_col(x, y=NULL, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[,0]
    ref <- crossprod(x0)
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x0, y=NULL, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x0, y=NULL, BPPARAM=MulticoreParam(2)))

    ref <- crossprod(x0, y)
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x0, y, BPPARAM=MulticoreParam(2)))

    y0 <- y[,0]
    ref <- crossprod(x, y0)
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x, y0, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    y0 <- y[0,]
    ref <- crossprod(x0, y0)
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_x_by_col(x0, y0, BPPARAM=MulticoreParam(2)))
})

test_that("bpcross_y_by_col works correctly", {
    x <- matrix(runif(1000), nrow=20)
    y <- matrix(runif(5000), nrow=20)
    ref <- crossprod(x, y)

    expect_equal(ref, BiocSingular:::bpcross_y_by_col(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross_y_by_col(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpcross_y_by_col(x, y, BPPARAM=MulticoreParam(3)))

    # Handles empties okay.
    x0 <- x[,0]
    ref <- crossprod(x0, y)
    expect_equivalent(ref, BiocSingular:::bpcross_y_by_col(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_y_by_col(x0, y, BPPARAM=MulticoreParam(2)))

    y0 <- y[,0]
    ref <- crossprod(x, y0)
    expect_equivalent(ref, BiocSingular:::bpcross_y_by_col(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_y_by_col(x, y0, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    y0 <- y[0,]
    ref <- crossprod(x0, y0)
    expect_equivalent(ref, BiocSingular:::bpcross_y_by_col(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_y_by_col(x0, y0, BPPARAM=MulticoreParam(2)))
})

test_that("bpcross_by_row works correctly", {
    x <- matrix(runif(1000), nrow=20)
    y <- matrix(runif(5000), nrow=20)
    ref <- crossprod(x, y)

    expect_equal(ref, BiocSingular:::bpcross_by_row(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross_by_row(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpcross_by_row(x, y, BPPARAM=MulticoreParam(3)))

    # Only one matrix.
    ref <- crossprod(x)
    expect_equal(ref, BiocSingular:::bpcross_by_row(x, y=NULL, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross_by_row(x, y=NULL, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpcross_by_row(x, y=NULL, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[,0]
    ref <- crossprod(x0)
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x0, y=NULL, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x0, y=NULL, BPPARAM=MulticoreParam(2)))

    ref <- crossprod(x0, y)
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x0, y, BPPARAM=MulticoreParam(2)))

    y0 <- y[,0]
    ref <- crossprod(x, y0)
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x, y0, BPPARAM=MulticoreParam(2)))

    x0 <- x[0,]
    y0 <- y[0,]
    ref <- crossprod(x0, y0)
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpcross_by_row(x0, y0, BPPARAM=MulticoreParam(2)))
})

test_that("bpcross overall function works correctly", {
    x <- matrix(runif(300), nrow=15)
    ref <- crossprod(x)
    expect_equal(ref, BiocSingular:::bpcross(x, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross(x, BPPARAM=MulticoreParam(2))) 
    expect_equal(ref, BiocSingular:::bpcross(x, BPPARAM=MulticoreParam(3))) 

    x <- matrix(runif(150), nrow=15)
    y <- matrix(runif(300), nrow=15)
    ref <- crossprod(x, y)
    expect_equal(ref, BiocSingular:::bpcross(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross(x, y, BPPARAM=MulticoreParam(2))) 
    expect_equal(ref, BiocSingular:::bpcross(x, y, BPPARAM=MulticoreParam(3))) 

    # Handles vector inputs.
    xv <- runif(20)
    ref <- crossprod(xv)
    expect_equal(ref, BiocSingular:::bpcross(xv, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross(xv, BPPARAM=MulticoreParam(2))) 

    y <- matrix(runif(300), nrow=20)
    ref <- crossprod(xv, y)
    expect_equal(ref, BiocSingular:::bpcross(xv, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross(xv, y, BPPARAM=MulticoreParam(2)))

    x <- matrix(runif(150), nrow=15)
    yv <- runif(15)
    ref <- crossprod(x, yv)
    expect_equal(ref, BiocSingular:::bpcross(x, yv, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpcross(x, yv, BPPARAM=MulticoreParam(2)))
})
