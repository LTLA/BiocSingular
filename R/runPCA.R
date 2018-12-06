#' @export
#' @importFrom DelayedArray DelayedArray sweep
#' @importFrom BiocGenerics nrow colMeans colSums t
runPCA <- function(x, rank, center=TRUE, scale=FALSE, BSPARAM=NULL, get.rotation=TRUE, get.pcs=TRUE, ...) {
    if (is.logical(center)) {
        if (center) {
            center <- colMeans(DelayedArray(x))
        } else {
            center <- NULL
        }
    }

    if (is.logical(scale)) {
        if (scale) {
            sub <- DelayedArray(x)
            if (!is.null(center)) {
                sub <- sweep(x, 2, center, "-")
            }
            scale <- colSums(sub^2) / (nrow(x) - 1L) # mimic scale() behaviour for any 'center'.
            scale <- sqrt(scale)
        } else {
            scale <- NULL
        }
    }

    svd.out <- runSVD(x, k=rank, 
        nu=ifelse(get.pcs, rank, 0),
        nv=ifelse(get.rotation, rank, 0),
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
