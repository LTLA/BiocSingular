#' @export
#' @importFrom utils head
#' @importFrom BiocParallel SerialParam
runExactSVD <- function(x, k=min(dim(x)), nu=k, nv=k, center=NULL, scale=NULL, fold=5L, BPPARAM=SerialParam())
# Wrapper for svd(), with options for faster calculation by taking the 
# cross-product for fat or tall matrices.
{
#    dx <- DelayedArray(x) # uncomment this line once *crossprod is supported by DAs.
    dx <- x
    if (!is.null(center)) {
        dx <- sweep(dx, 2, center, "-")
    }
    if (!is.null(scale)) {
        dx <- sweep(dx, 2, scale, "/")
    }

    if (nrow(x) > fold*ncol(x)) {
        y <- bpcross(dx, BPPARAM=BPPARAM)
        res <- .safe.svd(y, nu=0, nv=max(nu, nv))
        res$d <- sqrt(res$d)
        res$u <- sweep(dx %*% res$v[,seq_len(nu),drop=FALSE], 2, head(res$d, nu), "/")
        res$v <- res$v[,seq_len(nv),drop=FALSE]

    } else if (ncol(x) > fold*nrow(x)) {
        y <- bptcross(dx, BPPARAM=BPPARAM)
        res <- .safe.svd(y, nu=max(nu, nv), nv=0)
        res$d <- sqrt(res$d)
        v0 <- bpcross(dx, res$u[,seq_len(nv),drop=FALSE], BPPARAM=BPPARAM)
        res$v <- sweep(v0, 2, head(res$d, nv), "/")
        res$u <- res$u[,seq_len(nu),drop=FALSE]

    } else {
        res <- .safe.svd(as.matrix(dx), nu=nu, nv=nv)
    }

    res$d <- head(res$d, k)
    res
}

.safe.svd <- function(x, nu, nv) 
# Wrapper that guarantees return of U, D and V, even if nu or nv are zero. 
{
    out <- svd(x, nu=nu, nv=nv)
    if (!nu) {
        out$u <- matrix(0, nrow(x), 0)
    }
    if (!nv) {
        out$v <- matrix(0, ncol(x), 0)
    }
    out[c("d", "u", "v")]
}
