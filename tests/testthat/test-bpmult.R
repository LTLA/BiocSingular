# Tests bpmult() for parallelized matrix multiplication.
# library(testthat); library(BiocSingular); source("test-bpmult.R")

test_that("bpmult_x_by_row works correctly on two input matrices", {
    x <- matrix(runif(1000), ncol=20)
    y <- matrix(runif(5000), nrow=20)
    ref <- x%*%y

    expect_equal(ref, BiocSingular:::bpmult_x_by_row(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult_x_by_row(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult_x_by_row(x, y, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[0,]
    ref <- x0 %*% y
    expect_equal(ref, BiocSingular:::bpmult_x_by_row(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_row(x0, y, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_row(x0, y, BPPARAM=MulticoreParam(3)))

    y0 <- y[,0]
    ref <- x %*% y0
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_row(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_row(x, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_row(x, y0, BPPARAM=MulticoreParam(3)))

    # Both are empty.
    x0 <- x[,0]
    y0 <- y[0,]
    ref <- x0 %*% y0
    expect_equal(ref, BiocSingular:::bpmult_x_by_row(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_row(x0, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_row(x0, y0, BPPARAM=MulticoreParam(3)))
})

test_that("bpmult_x_by_col works correctly on two input matrices", {
    x <- matrix(runif(1000), ncol=25)
    y <- matrix(runif(5000), nrow=25)
    ref <- x%*%y

    expect_equal(ref, BiocSingular:::bpmult_x_by_col(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult_x_by_col(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult_x_by_col(x, y, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[0,]
    ref <- x0 %*% y
    expect_equal(ref, BiocSingular:::bpmult_x_by_col(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_col(x0, y, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_col(x0, y, BPPARAM=MulticoreParam(3)))

    y0 <- y[,0]
    ref <- x %*% y0
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_col(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_col(x, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_col(x, y0, BPPARAM=MulticoreParam(3)))

    # Both are empty.
    x0 <- x[,0]
    y0 <- y[0,]
    ref <- x0 %*% y0
    expect_equal(ref, BiocSingular:::bpmult_x_by_col(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_col(x0, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_x_by_col(x0, y0, BPPARAM=MulticoreParam(3)))
})

test_that("bpmult_y_by_col works correctly on two input matrices", {
    x <- matrix(runif(1000), ncol=50)
    y <- matrix(runif(5000), nrow=50)
    ref <- x%*%y

    expect_equal(ref, BiocSingular:::bpmult_y_by_col(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult_y_by_col(x, y, BPPARAM=MulticoreParam(2)))
    expect_equal(ref, BiocSingular:::bpmult_y_by_col(x, y, BPPARAM=MulticoreParam(3)))

    # One of these is empty.
    x0 <- x[0,]
    ref <- x0 %*% y
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x0, y, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x0, y, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x0, y, BPPARAM=MulticoreParam(3)))

    y0 <- y[,0]
    ref <- x %*% y0
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x, y0, BPPARAM=MulticoreParam(3)))

    # Both are empty.
    x0 <- x[,0]
    y0 <- y[0,]
    ref <- x0 %*% y0
    expect_equal(ref, BiocSingular:::bpmult_y_by_col(x0, y0, BPPARAM=SerialParam()))
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x0, y0, BPPARAM=MulticoreParam(2)))
    expect_equivalent(ref, BiocSingular:::bpmult_y_by_col(x0, y0, BPPARAM=MulticoreParam(3)))
})

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


test_that("bpcross overall function works correctly", {
    x <- matrix(runif(150), ncol=15)
    y <- matrix(runif(300), nrow=15)
    ref <- x %*% y
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=MulticoreParam(2))) # by row of 'x', which is divisible by 2.
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=MulticoreParam(3))) # by column, which is divisible by 3.
    expect_equal(ref, BiocSingular:::bpmult(x, y, BPPARAM=MulticoreParam(4))) # by row of 'y', which is divisible by 4.

    # Handles vector inputs.
    xv <- runif(20)
    y <- matrix(runif(300), nrow=20)
    ref <- xv %*% y
    expect_equal(ref, BiocSingular:::bpmult(xv, y, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(xv, y, BPPARAM=MulticoreParam(2)))

    x <- matrix(runif(150), nrow=15)
    yv <- runif(10)
    ref <- x %*% yv
    expect_equal(ref, BiocSingular:::bpmult(x, yv, BPPARAM=SerialParam()))
    expect_equal(ref, BiocSingular:::bpmult(x, yv, BPPARAM=MulticoreParam(2)))
})
