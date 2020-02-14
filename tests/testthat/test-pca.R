# Tests the runPCA() function.
# library(testthat); library(BiocSingular); source("setup.R"); source("test-pca.R")

set.seed(10000)
test_that("runPCA with exact SVD matches up to the reference implementation", {
    a <- matrix(rnorm(100000), ncol=20)
    out <- runPCA(a, rank=10)
    ref <- prcomp(a, rank.=10)

    expect_equal(out$sdev, head(ref$sdev, 10))
    expect_equal_besides_sign(out$rotation, ref$rotation)
    expect_equal_besides_sign(out$x, ref$x)

    # With scaling.
    out <- runPCA(a, rank=5, scale=TRUE)
    ref <- prcomp(a, rank.=5, scale=TRUE)

    expect_equal(out$sdev, head(ref$sdev, 5))
    expect_equal_besides_sign(out$rotation, ref$rotation)
    expect_equal_besides_sign(out$x, ref$x)

    # With scaling, but without centering.
    out <- runPCA(a, rank=5, scale=TRUE, center=FALSE)
    ref <- prcomp(a, rank.=5, scale=TRUE, center=FALSE)

    expect_equal(out$sdev, head(ref$sdev, 5))
    expect_equal_besides_sign(out$rotation, ref$rotation)
    expect_equal_besides_sign(out$x, ref$x)
})

set.seed(10001)
test_that("runPCA with approximate SVD (IRLBA) matches up to the reference implementation", {
    a <- matrix(rnorm(100000), ncol=50)
    set.seed(200)
    out <- runPCA(a, rank=10, BSPARAM=IrlbaParam(fold=Inf))
    set.seed(200)
    ref <- irlba::prcomp_irlba(a, n=10)

    expect_equal(out$sdev, ref$sdev)
    expect_equal_besides_sign(out$rotation, ref$rotation)
    expect_equal_besides_sign(out$x, ref$x)
})

set.seed(10002)
test_that("runPCA with randomized SVD matches up to the reference implementation", {
    a <- matrix(rnorm(100000), ncol=50)
    set.seed(200)
    out <- runPCA(a, rank=10, BSPARAM=RandomParam(fold=Inf))
    set.seed(200)
    ref <- rsvd::rpca(a, k=10, scale=FALSE)

    expect_equal(out$sdev, ref$sdev)
    expect_equal_besides_sign(out$rotation, ref$rotation)
    dimnames(out$x) <- NULL
    expect_equal_besides_sign(out$x, ref$x)
})

set.seed(500041)
test_that("runPCA handles named inputs", {
    y <- matrix(rnorm(10000), ncol=50)
    rownames(y) <- sprintf("THING_%i", seq_len(nrow(y)))
    colnames(y) <- sprintf("STUFF_%i", seq_len(ncol(y)))

    out <- runPCA(y, rank=5)
    expect_identical(rownames(out$x), rownames(y))
    expect_identical(rownames(out$rotation), colnames(y))
})
