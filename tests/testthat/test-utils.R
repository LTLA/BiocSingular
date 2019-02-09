# This explicitly unit tests utilities that are indirectly covered elsewhere.
# library(BiocSingular); library(testthat); source("test-utils.R")

test_that("use_crossprod works correctly", {
    expect_true(BiocSingular:::use_crossprod(matrix(0, nrow=10, ncol=100), 1))    
    expect_true(BiocSingular:::use_crossprod(matrix(0, nrow=10, ncol=100), 5))
    expect_true(BiocSingular:::use_crossprod(matrix(0, nrow=10, ncol=100), 10))
    expect_false(BiocSingular:::use_crossprod(matrix(0, nrow=10, ncol=100), 11))
    expect_false(BiocSingular:::use_crossprod(matrix(0, nrow=10, ncol=100), Inf))

    expect_true(BiocSingular:::use_crossprod(matrix(0, nrow=100, ncol=10), 1))
    expect_true(BiocSingular:::use_crossprod(matrix(0, nrow=100, ncol=10), 5))
    expect_true(BiocSingular:::use_crossprod(matrix(0, nrow=100, ncol=10), 10))
    expect_false(BiocSingular:::use_crossprod(matrix(0, nrow=100, ncol=10), 11))
    expect_false(BiocSingular:::use_crossprod(matrix(0, nrow=100, ncol=10), Inf))
})

test_that("standardize_matrix works correctly", {
    A <- matrix(runif(2000), 50, 40)

    # Cycling through all possible centering and scaling options.
    for (center in list(TRUE, FALSE, NULL, rnorm(ncol(A)))) {
        for (scale in list(TRUE, FALSE, NULL, 1+runif(ncol(A)))) {
            ref <- scale(A, 
                center=if (is.null(center)) FALSE else center, 
                scale=if (is.null(scale)) FALSE else scale)
            ref <- ref[,] # remove attributes

            out <- BiocSingular:::standardize_matrix(A, center=center, scale=scale)
            expect_s4_class(out, "DelayedArray")
            expect_equivalent(as.matrix(out), ref)
            
            out <- BiocSingular:::standardize_matrix(A, center=center, scale=scale, deferred=TRUE)
            expect_s4_class(out, "DeferredMatrix")
            expect_equivalent(as.matrix(out), ref)

            out <- BiocSingular:::standardize_matrix(A, center=center, scale=scale, deferred=TRUE, BPPARAM=BiocParallel::MulticoreParam(2))
            expect_s4_class(out, "DeferredMatrix")
            expect_s4_class(DelayedArray::seed(out)@.matrix, "DelayedArray")
            expect_equivalent(as.matrix(out), ref)
        }
    }
})

test_that("standardize_output_SVD works correctly", {
    out <- BiocSingular:::standardize_output_SVD(list(u=1:5, d=TRUE, v=1:10))
    expect_identical(names(out), c("d", "u", "v"))
    expect_identical(out$u, cbind(1:5))
    expect_identical(out$v, cbind(1:10))
    expect_identical(out$d, 1)

    # Handles alternative matrix types in 'u' and 'v'.
    U <- Matrix::rsparsematrix(10, 20, 0.12)
    V <- DelayedArray::DelayedArray(matrix(runif(200), 10, 20))
    out <- BiocSingular:::standardize_output_SVD(list(u=U, d=1:5, v=V))
    expect_identical(names(out), c("d", "u", "v"))
    expect_identical(out$u, as.matrix(U)) 
    expect_identical(out$v, as.matrix(V))
    expect_identical(out$d, as.numeric(1:5))
})
