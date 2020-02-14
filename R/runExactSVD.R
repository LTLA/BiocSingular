#' @export
#' @importFrom BiocParallel bpstart bpstop SerialParam
#' @importFrom utils head
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
runExactSVD <- function(x, k=min(dim(x)), nu=k, nv=k, center=FALSE, scale=FALSE, deferred=FALSE, fold=Inf, BPPARAM=SerialParam())
{
    checked <- check_numbers(x, k=k, nu=nu, nv=nv)
    k <- checked$k
    nv <- checked$nv
    nu <- checked$nu

    # Setting up the parallelization environment.
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    if (!.bpisup2(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred, BPPARAM=BPPARAM)
    if (use_crossprod(x, fold)) {
        res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=safe_svd)
    } else {
        res <- safe_svd(as.matrix(x), nu=nu, nv=nv)
        res$d <- head(res$d, k)
        res <- standardize_output_SVD(res, x)
    }

    res
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
    out
}
