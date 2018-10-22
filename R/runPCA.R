#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colMeans2 rowSums2
#' @importFrom BiocGenerics nrow
runPCA <- function(x, k, center=TRUE, scale=FALSE, BSPARAM=NULL, get.rotation=TRUE, get.pcs=TRUE, ...) {
    if (is.logical(center)) {
        if (center) {
            center <- colMeans2(DelayedArray(x))
        } else {
            center <- NULL
        }
    }

    if (is.logical(scale)) {
        if (scale) {
            scale <- rowSums2((t(DelayedArray(x)) - center)^2) / (nrow(x) - 1L) # mimic scale() behaviour for any 'center'.
            scale <- sqrt(scale)
        } else {
            scale <- NULL
        }
    }

    svd.out <- runSVD(x, k=k, 
        nu=ifelse(get.pcs, k, 0),
        nv=ifelse(get.rotation, k, 0),
        center=center, scale=scale, BSPARAM=BSPARAM, ...)

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
}
