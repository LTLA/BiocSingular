# Tests runExactSVD().
# library(testthat); library(BiocSingular); source("setup.R"); source("test-exact-svd.R")

set.seed(50000)
test_that("exact SVD works on square matrices", {
    y <- matrix(rnorm(10000), ncol=100, nrow=100)
    out <- runExactSVD(y, fold=Inf)
    ref <- svd(y)
    expect_equal(out, ref)

    # Handles truncation.
    out <- runExactSVD(y, k=10, nv=5, nu=3, fold=Inf)
    ref <- svd(y, nv=5, nu=3)
    ref$d <- ref$d[1:10]
    expect_equal(out, ref)

    out <- runExactSVD(y, k=1, nv=5, nu=3, fold=Inf)
    ref <- svd(y, nv=5, nu=3)
    ref$d <- ref$d[1]
    expect_equal(out, ref)

    # Handles zeroes.
    out0 <- runExactSVD(y, k=0, nv=0, nu=0)
    expect_equal(out0$d, numeric(0))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(out0$v, ref$v[,0])
})

set.seed(50001)
test_that("exact SVD works on thin matrices", {
    y <- matrix(rnorm(10000), ncol=10)
    out <- runExactSVD(y, fold=1)
    ref <- svd(y)
    expect_equal_svd(out, ref)

    # Handles truncation.
    out <- runExactSVD(y, k=5, nv=4, nu=2, fold=1)
    ref <- svd(y, nu=2, nv=4)
    ref$d <- ref$d[1:5]
    expect_equal_svd(out, ref)

    out <- runExactSVD(y, k=1, nv=6, nu=2, fold=1)
    ref <- svd(y, nu=2, nv=6)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)

    # Handles zeroes.
    out0 <- runExactSVD(y, k=0, nv=0, nu=0, fold=1)
    expect_equal(out0$d, numeric(0))
    expect_equal(out0$u, ref$u[,0])
    expect_equal(out0$v, ref$v[,0])
})

set.seed(50002)
test_that("exact SVD works on fat matrices", {
    y <- matrix(rnorm(10000), nrow=10)
    out <- runExactSVD(y, fold=1)
    ref <- svd(y)
    expect_equal_svd(out, ref)

    # Handles truncation.
    out <- runExactSVD(y, k=7, nv=6, nu=2, fold=1)
    ref <- svd(y, nu=2, nv=6)
    ref$d <- ref$d[1:7]
    expect_equal_svd(out, ref)

    out <- runExactSVD(y, k=1, nv=6, nu=2, fold=1)
    ref <- svd(y, nu=2, nv=6)
    ref$d <- ref$d[1]
    expect_equal_svd(out, ref)

    # Handles zeroes.
    out0 <- runExactSVD(y, k=0, nv=0, nu=0, fold=1)
    expect_equal(out0$d, numeric(0))
    expect_equivalent(out0$u, ref$u[,0])
    expect_equivalent(out0$v, ref$v[,0])
})

set.seed(50003)
test_that("exact SVD works with centering and scaling", {
    y <- matrix(rnorm(10000), ncol=50)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))

    ry <- scale(y, center=center, scale=scale)
    ref <- svd(ry)
    out <- runExactSVD(y, center=center, scale=scale, fold=Inf)
    expect_equal_svd(out, ref)

    # Works with the crossproduct.
    y <- matrix(rnorm(10000), ncol=10)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))

    ry <- scale(y, center=center, scale=scale)
    ref <- svd(ry)
    out <- runExactSVD(y, center=center, scale=scale, fold=1)
    expect_equal_svd(out, ref)

    # Works with deferred operations.
    y <- matrix(rnorm(10000), ncol=10)
    center <- runif(ncol(y))
    scale <- runif(ncol(y))

    ry <- scale(y, center=center, scale=scale)
    ref <- svd(ry)
    out <- runExactSVD(y, center=center, scale=scale, deferred=TRUE)
    expect_equal_svd(out, ref)
})

set.seed(500041)
test_that("exact SVD handles named inputs", {
    y <- matrix(rnorm(10000), ncol=50)
    rownames(y) <- sprintf("THING_%i", seq_len(nrow(y)))
    colnames(y) <- sprintf("STUFF_%i", seq_len(ncol(y)))

    out <- runExactSVD(y, fold=Inf)
    expect_identical(rownames(out$u), rownames(y))
    expect_identical(rownames(out$v), colnames(y))

    out <- runExactSVD(y, fold=1)
    expect_identical(rownames(out$u), rownames(y))
    expect_identical(rownames(out$v), colnames(y))
})

set.seed(50004)
test_that("exact SVD fails gracefully with silly inputs", {
    y <- matrix(rnorm(10000), ncol=50)
    expect_error(runExactSVD(y, k=-1), "non-negative")
    expect_error(runExactSVD(y, nu=-1), "non-negative")
    expect_error(runExactSVD(y, nv=-1), "non-negative")

    expect_warning(runExactSVD(y, k=1e6), "requested than available")
    expect_warning(runExactSVD(y, nu=1e6), "requested than available")
    expect_warning(runExactSVD(y, nv=1e6), "requested than available")

    expect_warning(out <- runExactSVD(y, k=Inf), NA)
    expect_identical(length(out$d), min(dim(y)))
    expect_warning(runExactSVD(y, nu=Inf), NA)
    expect_identical(ncol(out$u), min(dim(y)))
    expect_warning(runExactSVD(y, nv=Inf), NA)
    expect_identical(ncol(out$v), min(dim(y)))
})
