# This explicitly unit tests utilities that are indirectly covered elsewhere.
# library(BiocSingular); library(testthat); source("setup.R"); source("test-utils.R")

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
    for (n in 1:2) {
        helper <- function(...) {
            BiocSingular:::.compute_center_and_scale(..., scale=TRUE, nthreads=n)$scale
        }

        for (it in 1:4) {
            if (it==1L) {
                A <- matrix(runif(2000), 50, 40)
            } else if (it==2L) {
                A <- t(DelayedArray(matrix(rpois(2000, lambda=5), 50, 40))) # dense, prefers rows
            } else if (it==3L) {
                A <- Matrix::rsparsematrix(50, 40, density=0.1)
            } else {
                A <- t(DelayedArray(Matrix::rsparsematrix(50, 40, density=0.1))) # sparse, prefers rows
            }

            out <- helper(A, center=FALSE)
            A0 <- as.matrix(A)
            ref <- sqrt(Matrix::colSums(A0^2)/(nrow(A0)-1))
            expect_equal(out, ref)

            center <- rnorm(ncol(A0))
            out <- helper(A, center=center)
            B <- sweep(A0, 2, center, "-")
            ref <- sqrt(Matrix::colSums(B^2)/(nrow(B)-1))
            expect_equal(out, ref)

            out <- helper(A, center=TRUE)
            center <- colMeans(A0)
            B <- sweep(A0, 2, center, "-")
            ref <- sqrt(Matrix::colSums(B^2)/(nrow(B)-1))
            expect_equal(out, ref)

            expect_identical(helper(A[,0], numeric(0)), numeric(0))
            expect_identical(helper(A[0,], TRUE), rep(NA_real_, ncol(A)))
            expect_identical(helper(A[1,,drop=FALSE], TRUE), rep(NA_real_, ncol(A)))
            expect_identical(helper(A[0,], NULL), rep(NA_real_, ncol(A)))
            expect_error(helper(A, 1), "should be equal")
        }
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
            expect_s4_class(out, "ScaledMatrix")
            expect_equivalent(as.matrix(out), ref)

            if (.Platform$OS.type!="windows") { # as safeBPParam uses SerialParam, which doesn't wrap the matrix in a DA.
                out <- BiocSingular:::standardize_matrix(A, center=center, scale=scale, deferred=TRUE, BPPARAM=safeBPParam(2))
                expect_s4_class(out, "ScaledMatrix")
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
