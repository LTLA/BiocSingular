# Tests runRandomSVD().
# library(testthat); library(BiocSingular); source("setup.R"); source("test-random-svd.R")

library(rsvd)

set.seed(80000)
test_that("Random SVD works on input matrices", {
    y <- matrix(rnorm(50000), ncol=200)
    set.seed(100)
    out <- runRandomSVD(y, k=5, fold=Inf)
    set.seed(100)
    ref <- rsvd(y, k=5, nu=5, nv=5)
    expect_equal_svd(out, ref[c("d", "u", "v")])

    # Handles truncation.
    set.seed(100)
    out <- runRandomSVD(y, k=10, nv=5, nu=3, fold=Inf)
    set.seed(100)
    ref <- rsvd(y, k=10, nv=5, nu=3)
    ref$v <- ref$v[,1:5]
    expect_equal_svd(out, ref[c("d", "u", "v")])

    set.seed(100)
    out <- runRandomSVD(y, k=1, nv=5, nu=3, fold=Inf)
    set.seed(100)
    ref <- rsvd(y, k=5, nv=5, nu=3)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref[c("d", "u", "v")])
})

set.seed(8000001)
test_that("Random SVD handles zeroes", {
    y <- matrix(rnorm(50000), ncol=200)
    ref <- rsvd(y, k=5, nu=5, nv=5)

    out0 <- runRandomSVD(y, k=0, nv=0, nu=0)
    expect_equal(out0$d, numeric(0))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(out0$v, ref$v[,0])

    out0 <- runRandomSVD(y, k=5, nv=5, nu=0)
    expect_equal(length(out0$d), length(ref$d))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(dim(out0$v), dim(ref$v))

    out0 <- runRandomSVD(y, k=5, nv=0, nu=5)
    expect_equal(length(out0$d), length(ref$d))
    expect_equal(dim(out0$u), dim(ref$u))
    expect_equal(out0$v, ref$v[,0])
})

set.seed(80001)
test_that("Random SVD works on thin matrices", {
    y <- matrix(rnorm(10000), ncol=10)
    set.seed(200)
    out <- runRandomSVD(y, k=3, fold=1)
    set.seed(200)
    ref <- rsvd(y, k=3, nv=3, nu=3)
    expect_equal_svd(out, ref)

    # Handles truncation.
    set.seed(200)
    out <- runRandomSVD(y, k=5, nv=3, nu=2, fold=1)
    set.seed(200)
    ref <- rsvd(y, k=5, nu=2, nv=3)
    ref$v <- ref$v[,1:3]
    expect_equal_svd(out, ref)

    set.seed(200)
    out <- runRandomSVD(y, k=1, nv=3, nu=2, fold=1)
    set.seed(200)
    ref <- rsvd(y, k=3, nu=2, nv=3)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)
})

set.seed(80002)
test_that("Random SVD works on fat matrices", {
    y <- matrix(rnorm(10000), nrow=10)
    set.seed(300)
    out <- runRandomSVD(y, k=4, fold=1)
    set.seed(300)
    ref <- rsvd(y, k=4, nu=4, nv=4)
    expect_equal_svd(out, ref)

    # Handles truncation.
    set.seed(300)
    out <- runRandomSVD(y, k=5, nv=4, nu=2, fold=1)
    set.seed(300)
    ref <- rsvd(y, k=5, nu=2, nv=4)
    expect_equal_svd(out, ref)

    set.seed(300)
    out <- runRandomSVD(y, k=1, nv=6, nu=2, fold=1)
    set.seed(300)
    ref <- rsvd(y, k=6, nu=2, nv=6)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)
})

set.seed(80003)
test_that("Random SVD works with parallelization", {
    y <- matrix(rnorm(50000), ncol=100)
    set.seed(100)
    ref <- runRandomSVD(y, k=5, fold=Inf, p=50, q=20)
    set.seed(100)
    out <- runRandomSVD(y, k=5, BPPARAM=safeBPParam(2), fold=Inf, p=50, q=20)
    expect_equal_svd(ref, out, tol=1e-6)

    # Handles truncation.
    set.seed(100)
    ref <- runRandomSVD(y, k=5, nv=2, nu=3, fold=Inf, p=50, q=20)
    set.seed(100)
    out <- runRandomSVD(y, k=5, nv=2, nu=3, BPPARAM=safeBPParam(2), fold=Inf, p=50, q=20)
    expect_equal_svd(out, ref[c("d", "u", "v")], tol=1e-6)
})

set.seed(80004)
test_that("Random SVD works with centering and scaling", {
    y <- matrix(rnorm(10000), ncol=50)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))

    set.seed(100)
    out <- runRandomSVD(y, k=5, center=center, scale=scale, fold=Inf)
    set.seed(100)
    ry <- scale(y, center=center, scale=scale)
    ref <- rsvd(ry, k=5, nv=5, nu=5)
    expect_equal_svd(out, ref)

    # Works with the cross-product.
    y <- matrix(rnorm(10000), ncol=10)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))
   
    ry <- scale(y, center=center, scale=scale)
    set.seed(200)
    ref <- rsvd(ry, k=6, nu=6, nv=6)
    set.seed(200)
    out <- runRandomSVD(y, k=6, center=center, scale=scale, fold=1)
    expect_equal_svd(out, ref)

    # Works with the deferred operations. 
    set.seed(100)
    ref <- runRandomSVD(ry, k=7)
    set.seed(100)
    out <- runRandomSVD(y, k=7, center=center, scale=scale, deferred=TRUE)
    expect_equal_svd(out, ref)
})

set.seed(800041)
test_that("Random SVD handles named inputs", {
    y <- matrix(rnorm(10000), ncol=50)
    rownames(y) <- sprintf("THING_%i", seq_len(nrow(y)))
    colnames(y) <- sprintf("STUFF_%i", seq_len(ncol(y)))

    out <- runRandomSVD(y, k=3, fold=Inf)
    expect_identical(rownames(out$u), rownames(y))
    expect_identical(rownames(out$v), colnames(y))

    out <- runRandomSVD(y, k=3, fold=1)
    expect_identical(rownames(out$u), rownames(y))
    expect_identical(rownames(out$v), colnames(y))
})

set.seed(80005)
test_that("Random SVD fails gracefully with silly inputs", {
    y <- matrix(rnorm(10000), ncol=50)
    expect_error(runRandomSVD(y, k=-1), "non-negative")
    expect_error(runRandomSVD(y, nu=-1), "non-negative")
    expect_error(runRandomSVD(y, nv=-1), "non-negative")

    expect_warning(runRandomSVD(y, k=1e6), "requested than available")
    expect_warning(runRandomSVD(y, nu=1e6), "requested than available")
    expect_warning(runRandomSVD(y, nv=1e6), "requested than available")
})
