#' @export
setMethod("runPCA", "ANY", function(x, rank, center=TRUE, scale=FALSE, get.rotation=TRUE, get.pcs=TRUE, ...) 
# Converts SVD results to PCA results.
{
    svd.out <- runSVD(x, k=rank, 
        nu=ifelse(get.pcs, rank, 0),
        nv=ifelse(get.rotation, rank, 0),
        center=center, scale=scale, ...)

    # Not naming sdev for consistency with prcomp().
    out <- list(sdev=svd.out$d / sqrt(nrow(x) - 1))

    NAMEFUN <- function(n) sprintf("PC%i", seq_len(n))

    if (get.rotation) {
        out$rotation <- svd.out$v
        colnames(out$rotation) <- NAMEFUN(ncol(out$rotation))
    } 
    if (get.pcs) {
        out$x <- sweep(svd.out$u, 2, svd.out$d, "*")
        colnames(out$x) <- NAMEFUN(ncol(out$x))
    }

    out
})

#' @export
setMethod("runPCA", "ResidualMatrix", function(x, rank, center=TRUE, scale=FALSE, get.rotation=TRUE, get.pcs=TRUE, ...) {
    if (center && is_centered(seed(x))) {
        center <- FALSE
    }
    callNextMethod()
})
