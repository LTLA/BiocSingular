# Tests the bs_matrix implementation.
# library(testthat); library(BiocSingular); source("test-bs-matrix.R")

##########################

test_that("bs_matrix right multiplication works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(100000), ncol=20)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)

        # Multiply by a vector.
        z <- rnorm(ncol(y))
        expect_equal(bs.y %*% z, ref.y %*% z)

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(y)*10), ncol=10)
        expect_equal(bs.y %*% z, ref.y %*% z)

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=ncol(y))
        expect_equal(bs.y %*% z, ref.y %*% z)
    }
})

test_that("bs_matrix left multiplication works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(100000), ncol=20)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)

        # Multiply by a vector.
        z <- rnorm(nrow(y))
        expect_equal(z %*% bs.y, z %*% ref.y)

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(y)*10), nrow=10)
        expect_equal(z %*% bs.y, z %*% ref.y)

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=nrow(y))
        expect_equal(z %*% bs.y, z %*% ref.y)
    }
})

test_that("bs_matrix dual multiplication fails as expected", {
    y <- matrix(rnorm(400), ncol=20)
    bs.y <- BiocSingular:::bs_matrix(y, NULL, NULL)
    expect_error(bs.y %*% bs.y, "not yet supported")
})

##########################

test_that("bs_matrix lonely crossproduct works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(10000), ncol=200)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)
        expect_equal(BiocSingular:::crossprod(bs.y), crossprod(ref.y))
    }
})

test_that("bs_matrix left crossproduct works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(10000), ncol=200)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)

        # Multiply by a vector.
        z <- rnorm(nrow(y))
        expect_equal(BiocSingular:::crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(y)*10), ncol=10)
        expect_equal(BiocSingular:::crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(y))
        expect_equal(BiocSingular:::crossprod(bs.y, z), crossprod(ref.y, z))
    }
})

test_that("bs_matrix right crossproduct works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(10000), ncol=200)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)

        # Multiply by a vector.
        z <- rnorm(nrow(y))
        expect_equal(BiocSingular:::crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(y)*10), ncol=10)
        expect_equal(BiocSingular:::crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(y))
        expect_equal(BiocSingular:::crossprod(z, bs.y), crossprod(z, ref.y))
    }
})

##########################

test_that("bs_matrix lonely tcrossproduct works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(10000), ncol=200)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)
        expect_equal(BiocSingular:::tcrossprod(bs.y), tcrossprod(ref.y))
    }
})

test_that("bs_matrix left tcrossproduct works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(10000), ncol=200)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)

        # Multiply by a vector (this doesn't work).
        z <- rnorm(ncol(y))
        expect_error(BiocSingular:::tcrossprod(bs.y, z), "non-conformable")
        expect_error(tcrossprod(ref.y, z), "non-conformable")

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(y)*10), nrow=10)
        expect_equal(BiocSingular:::tcrossprod(bs.y, z), tcrossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(y))
        expect_equal(BiocSingular:::tcrossprod(bs.y, z), tcrossprod(ref.y, z))
    }
})

test_that("bs_matrix right crossproduct works as expected", {
    for (it in 1:4) {
        y <- matrix(rnorm(10000), ncol=200)
        center <- scale <- NULL

        if (it==1L) {
            center <- colMeans(y)
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=center, scale=scale)
        } else if (it==2L) {
            center <- rnorm(ncol(y))
            ref.y <- scale(y, center=center, scale=FALSE)
        } else if (it==3L) {
            scale <- runif(ncol(y))
            ref.y <- scale(y, center=FALSE, scale=scale)
        } else {
            ref.y <- y
        }

        bs.y <- BiocSingular:::bs_matrix(y, center, scale)

        # Multiply by a vector.
        z <- rnorm(ncol(y))
        expect_equal(BiocSingular:::tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(y)*10), nrow=10)
        expect_equal(BiocSingular:::tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(y))
        expect_equal(BiocSingular:::tcrossprod(z, bs.y), tcrossprod(z, ref.y))
    }
})
