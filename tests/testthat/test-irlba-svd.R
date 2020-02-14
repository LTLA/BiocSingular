# Tests runIrlbaSVD().
# library(testthat); library(BiocSingular); source("setup.R"); source("test-irlba-svd.R")

library(irlba)

set.seed(9000)
test_that("IRLBA works on input matrices", {

    y <- matrix(rnorm(50000), ncol=200)
    set.seed(100)
    out <- runIrlbaSVD(y, k=5, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5)
    expect_equal_svd(out, ref[c("d", "u", "v")])

    # Handles truncation.
    set.seed(100)
    out <- runIrlbaSVD(y, k=10, nv=5, nu=3, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=10, nu=3)
    ref$v <- ref$v[,1:5]
    expect_equal_svd(out, ref[c("d", "u", "v")])

    set.seed(100)
    out <- runIrlbaSVD(y, k=1, nv=5, nu=3, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5, nu=3)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref[c("d", "u", "v")])

    # Handles zeroes.
    out0 <- runIrlbaSVD(y, k=0, nv=0, nu=0)
    expect_equal(out0$d, numeric(0))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(out0$v, ref$v[,0])
})

set.seed(9001)
test_that("IRLBA works on thin matrices", {
    y <- matrix(rnorm(10000), ncol=10)
    set.seed(200)
    out <- runIrlbaSVD(y, k=3, fold=1)
    set.seed(200)
    ref <- irlba(y, nv=3, nu=3)
    expect_equal_svd(out, ref)

    # Handles truncation.
    set.seed(200)
    out <- runIrlbaSVD(y, k=5, nv=3, nu=2, fold=1)
    set.seed(200)
    ref <- irlba(y, nu=2, nv=5)
    ref$v <- ref$v[,1:3]
    expect_equal_svd(out, ref)

    set.seed(200)
    out <- runIrlbaSVD(y, k=1, nv=3, nu=2, fold=1)
    set.seed(200)
    ref <- irlba(y, nu=2, nv=3)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)
})

set.seed(9002)
test_that("IRLBA works on fat matrices", {
    y <- matrix(rnorm(10000), nrow=10)
    set.seed(300)
    out <- runIrlbaSVD(y, k=4, fold=1)
    set.seed(300)
    ref <- irlba(y, nu=4, nv=4)
    expect_equal_svd(out, ref)

    # Handles truncation.
    set.seed(300)
    out <- runIrlbaSVD(y, k=5, nv=4, nu=2, fold=1)
    set.seed(300)
    ref <- irlba(y, nu=2, nv=5)
    ref$v <- ref$v[,1:4]
    expect_equal_svd(out, ref)

    set.seed(300)
    out <- runIrlbaSVD(y, k=1, nv=6, nu=2, fold=1)
    set.seed(300)
    ref <- irlba(y, nu=2, nv=6)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)
})

set.seed(9003)
test_that("IRLBA works with alternative multiplication", {
    # Lower tol to avoid numeric precision issues.
    y <- matrix(rnorm(50000), ncol=100)
    set.seed(100)
    out <- runIrlbaSVD(y, k=5, BPPARAM=safeBPParam(2), tol=1e-8, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5, tol=1e-8)
    expect_equal_svd(out, ref[c("d", "u", "v")], tol=1e-6)

    # Handles truncation.
    set.seed(100)
    out <- runIrlbaSVD(y, k=5, nv=2, nu=3, BPPARAM=safeBPParam(2), tol=1e-8, fold=Inf)
    set.seed(100)
    ref <- irlba(y, nv=5, nu=3, tol=1e-8)
    ref$v <- ref$v[,1:2]
    expect_equal_svd(out, ref[c("d", "u", "v")], tol=1e-6)
})

set.seed(9004)
test_that("IRLBA works with centering and scaling", {
    y <- matrix(rnorm(10000), ncol=50)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))

    set.seed(100)
    ref <- irlba(y, k=5, center=center, scale=scale, tol=1e-8)
    set.seed(100)
    out <- runIrlbaSVD(y, k=5, center=center, scale=scale, fold=Inf, tol=1e-8) 
    expect_equal_svd(out, ref)

    ry <- scale(y, center=center, scale=scale)
    set.seed(100)
    ref2 <- irlba(ry, nv=5, nu=5, tol=1e-8)
    expect_equal_svd(out, ref2)

    # Works with logical values.
    ry <- scale(y, center=TRUE, scale=TRUE)
    set.seed(300)
    ref <- irlba(ry, nu=5, nv=5, tol=1e-8)
    set.seed(300)
    out <- runIrlbaSVD(y, k=5, center=TRUE, scale=TRUE, fold=1, tol=1e-8)
    expect_equal_svd(out, ref, tol=1e-6)
})

set.seed(9004)
test_that("IRLBA centering and scaling interact happily with other modes", {
    # Works with the cross-product.
    y <- matrix(rnorm(10000), ncol=20)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))
    
    ry <- scale(y, center=center, scale=scale)
    set.seed(200)
    ref <- irlba(ry, nu=5, nv=5, tol=1e-8)
    set.seed(200)
    out <- runIrlbaSVD(y, k=5, center=center, scale=scale, fold=1, tol=1e-8)
    expect_equal_svd(out, ref, tol=1e-6)

    # Works with the alternative multiplication.
    set.seed(100)
    out <- runIrlbaSVD(y, k=6, center=center, scale=scale, BPPARAM=safeBPParam(2), tol=1e-8, fold=Inf)
    set.seed(100)
    ref <- irlba(ry, nu=6, nv=6, tol=1e-8)
    expect_equal_svd(out, ref, tol=1e-6)

    # Works with our deferred multiplication (which also requires parallelization,
    # in order to not just use irlba's deferred methods directly).
    set.seed(100)
    ref <- irlba(ry, nu=7, nv=7, tol=1e-8)
    set.seed(100)
    out <- runIrlbaSVD(y, k=7, center=center, scale=scale, deferred=TRUE, BPPARAM=safeBPParam(2), tol=1e-8, fold=Inf)
    expect_equal_svd(out, ref, tol=1e-6)
})

set.seed(90041)
test_that("irlba SVD handles named inputs", {
    y <- matrix(rnorm(10000), ncol=50)
    rownames(y) <- sprintf("THING_%i", seq_len(nrow(y)))
    colnames(y) <- sprintf("STUFF_%i", seq_len(ncol(y)))

    out <- runIrlbaSVD(y, k=3, fold=Inf)
    expect_identical(rownames(out$u), rownames(y))
    expect_identical(rownames(out$v), colnames(y))

    out <- runIrlbaSVD(y, k=3, fold=1)
    expect_identical(rownames(out$u), rownames(y))
    expect_identical(rownames(out$v), colnames(y))
})

set.seed(9005)
test_that("IRLBA fails gracefully with silly inputs", {
    y <- matrix(rnorm(10000), ncol=50)
    expect_error(runIrlbaSVD(y, k=-1), "non-negative")
    expect_error(runIrlbaSVD(y, nu=-1), "non-negative")
    expect_error(runIrlbaSVD(y, nv=-1), "non-negative")

    expect_warning(runIrlbaSVD(y, k=1e6), "requested than available")
    expect_warning(runIrlbaSVD(y, nu=1e6), "requested than available")
    expect_warning(runIrlbaSVD(y, nv=1e6), "requested than available")
})
