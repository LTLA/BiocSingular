#' @importFrom BiocParallel SerialParam
#' @importFrom utils head
#' @importFrom BiocGenerics nrow ncol
svd_via_crossprod <- function(x, k, nu=k, nv=k, FUN=svd, ..., BPPARAM=SerialParam()) 
# Computes the SVD via a crossproduct, using 'FUN' with arguments '...'. 
# We assume that any centering/scaling of 'x' has already been applied to 'x'.
# For this reason, we can't easily run 'FUN(x)' directly, as this might require centering/scaling within 'FUN'.
{
    if (nrow(x) > ncol(x)) {
        y <- bpcross(x, BPPARAM=BPPARAM)
        res <- FUN(y, nu=0, nv=max(nu, nv, k), ...)
        res$d <- sqrt(res$d)

        u0 <- bpmult(x, res$v[,seq_len(nu),drop=FALSE], BPPARAM=BPPARAM)
        res$u <- sweep(u0, 2, head(res$d, nu), "/")
        res$v <- res$v[,seq_len(nv),drop=FALSE]

    } else {
        y <- bptcross(x, BPPARAM=BPPARAM)
        res <- FUN(y, nu=max(nu, nv, k), nv=0, ...)
        res$d <- sqrt(res$d)

        v0 <- bpcross(x, res$u[,seq_len(nv),drop=FALSE], BPPARAM=BPPARAM)
        res$v <- sweep(v0, 2, head(res$d, nv), "/")
        res$u <- res$u[,seq_len(nu),drop=FALSE]
    }

    res$d <- head(res$d, k)
	return(res)
}

#' @importFrom BiocGenerics nrow ncol
use_crossprod <- function(x, fold) {
    if (any(dim(x)==0L)) { # avoid problems when 'fold=Inf'.
        return(FALSE)
    }
    return(nrow(x) >= fold*ncol(x) || ncol(x) >= nrow(x)*fold)
}
