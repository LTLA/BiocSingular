#' @export
#' @importFrom BiocParallel bpstart bpstop bpisup bpparam register
#' @importFrom utils head
runExactSVD <- function(x, k=min(dim(x)), nu=k, nv=k, center=NULL, scale=NULL, deferred=FALSE, fold=5L, BPPARAM=NULL)
# Wrapper for svd(), with options for faster calculation by taking the 
# cross-product for fat or tall matrices.
{
    # Setting up the parallelization environment.
    if (is.null(BPPARAM)) {
        BPPARAM <- bpparam()
    } else {
        old <- bpparam()
        register(BPPARAM)
        on.exit(register(old))
    }
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred)
    if (use_crossprod(x, fold)) {
        res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=safe_svd)
    } else {
        res <- safe_svd(as.matrix(x), nu=nu, nv=nv)
        res$d <- head(res$d, k)
    }

    return(res)
}

safe_svd <- function(x, nu, nv) 
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
