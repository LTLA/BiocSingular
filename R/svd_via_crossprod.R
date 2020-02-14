#' @importFrom BiocParallel SerialParam
#' @importFrom utils head
#' @importFrom BiocGenerics nrow ncol
#' @importFrom DelayedArray sweep
svd_via_crossprod <- function(x, k, nu=k, nv=k, FUN=svd, ...) 
# Computes the SVD via a crossproduct, using 'FUN' with arguments '...'. 
# We assume that any centering/scaling has already been applied to 'x'.
{
    if (nrow(x) > ncol(x)) {
        y <- as.matrix(crossprod(x))
        res <- FUN(y, nu=0, nv=max(nu, nv, k), ...)
        res$d <- sqrt(res$d)

        u0 <- x %*% res$v[,seq_len(nu),drop=FALSE]
        res$u <- sweep(u0, 2, head(res$d, nu), "/")
        res$v <- as.matrix(res$v[,seq_len(nv),drop=FALSE])

    } else {
        y <- as.matrix(tcrossprod(x))
        res <- FUN(y, nu=max(nu, nv, k), nv=0, ...)
        res$d <- sqrt(res$d)

        v0 <- crossprod(x, res$u[,seq_len(nv),drop=FALSE])
        res$v <- sweep(v0, 2, head(res$d, nv), "/")
        res$u <- as.matrix(res$u[,seq_len(nu),drop=FALSE])
    }

    res$d <- head(res$d, k)
    standardize_output_SVD(res, x)
}

#' @importFrom BiocGenerics nrow ncol
use_crossprod <- function(x, fold) {
    if (any(dim(x)==0L)) { # avoid problems when 'fold=Inf'.
        FALSE
    } else {
        nrow(x) >= fold*ncol(x) || ncol(x) >= nrow(x)*fold
    }
}
