expect_equal_besides_sign <- function(left, right, ...) {
    ratio <- left/right
    right2 <- sweep(right, 2, sign(ratio[1,]), FUN="*")
    expect_equal(left, right2, ...)
}

expect_equal_svd <- function(left, right, ...) {
    expect_equal_besides_sign(left$u, right$u, ...)
    expect_equal_besides_sign(left$v, right$v, ...)
    expect_equal(left$d, right$d, ...)
}

expect_equal_product <- function(x, y) {
    expect_s4_class(x, "DelayedMatrix")
    X <- as.matrix(x)

    # standardize NULL dimnames.
    if (all(lengths(dimnames(X))==0L)) dimnames(X) <- NULL 
    if (all(lengths(dimnames(y))==0L)) dimnames(y) <- NULL 
    expect_equal(X, y)
}

# Because SnowParam() is too slow, yet MulticoreParam() fails on Windows.
# See discussion at https://github.com/Bioconductor/BiocParallel/issues/98.
safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
