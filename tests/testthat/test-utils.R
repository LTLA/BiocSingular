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

test_that("scale calculations work correctly", {
    for (it in 1:4) {
        if (it==1L) {
            A <- matrix(runif(2000), 50, 40)
        } else if (it==2L) {
            A <- matrix(rpois(2000, lambda=5), 50, 40)
        } else if (it==3L) {
            A <- Matrix::rsparsematrix(50, 40, density=0.1)
        } else {
            A <- as(Matrix::rsparsematrix(50, 40, density=0.1), "dgTMatrix")
        }

        out <- BiocSingular:::compute_scale(A, NULL)
        ref <- sqrt(Matrix::colSums(A^2)/(nrow(A)-1))
        expect_equal(out, ref)

        center <- rnorm(ncol(A))
        out <- BiocSingular:::compute_scale(A, center)
        B <- sweep(A, 2, center, "-")
        ref <- sqrt(Matrix::colSums(B^2)/(nrow(B)-1))
        expect_equal(out, ref)

        expect_identical(BiocSingular:::compute_scale(A[,0], numeric(0)), numeric(0))
        expect_identical(BiocSingular:::compute_scale(A[0,], NULL), rep(NA_real_, ncol(A)))
        expect_error(BiocSingular:::compute_scale(A, 1), "should be equal")
    }
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

            if (.Platform$OS.type!="windows") { # as safeBPParam uses SerialParam, which doesn't wrap the matrix in a DA.
                out <- BiocSingular:::standardize_matrix(A, center=center, scale=scale, deferred=TRUE, BPPARAM=safeBPParam(2))
                expect_s4_class(out, "DeferredMatrix")
                expect_s4_class(DelayedArray::seed(out)@.matrix, "DelayedArray")
                expect_equivalent(as.matrix(out), ref)
            }
        }
    }
})

test_that("standardize_output_SVD works correctly", {
    out <- BiocSingular:::standardize_output_SVD(list(u=1:5, d=TRUE, v=1:10), NULL)
    expect_identical(names(out), c("d", "u", "v"))
    expect_identical(out$u, cbind(1:5))
    expect_identical(out$v, cbind(1:10))
    expect_identical(out$d, 1)

    # Handles alternative matrix types in 'u' and 'v'.
    U <- Matrix::rsparsematrix(10, 20, 0.12)
    V <- DelayedArray::DelayedArray(matrix(runif(200), 10, 20))
    out <- BiocSingular:::standardize_output_SVD(list(u=U, d=1:5, v=V), NULL)
    expect_identical(names(out), c("d", "u", "v"))
    expect_identical(out$u, as.matrix(U)) 
    expect_identical(out$v, as.matrix(V))
    expect_identical(out$d, as.numeric(1:5))

    # Sets names, if the matrix was named.
    mat <- matrix(0, 5, 10, dimnames=list(letters[1:5], LETTERS[1:10]))
    out <- BiocSingular:::standardize_output_SVD(list(u=1:5, d=TRUE, v=1:10), mat)
    expect_identical(rownames(out$u), rownames(mat))
    expect_identical(rownames(out$v), colnames(mat))
})
