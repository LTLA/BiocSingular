#' @export
setMethod("runPCA", "ANY", function(x, rank, center=TRUE, scale=FALSE, get.rotation=TRUE, get.pcs=TRUE, ...) 
# Converts SVD results to PCA results.
{
    svd.out <- runSVD(x, k=rank, 
        nu=ifelse(get.pcs, rank, 0),
        nv=ifelse(get.rotation, rank, 0),
        center=center, scale=scale, ...)

    out <- list(sdev=svd.out$d / sqrt(nrow(x) - 1))
    if (get.rotation) {
        out$rotation <- svd.out$v
        colnames(out$rotation) <- sprintf("PC%i", seq_len(ncol(out$rotation)))
    } 
    if (get.pcs) {
        out$x <- sweep(svd.out$u, 2, svd.out$d, "*")
        colnames(out$x) <- sprintf("PC%i", seq_len(ncol(out$x)))
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
