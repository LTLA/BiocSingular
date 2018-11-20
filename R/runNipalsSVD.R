#' @export
#' @importFrom BiocParallel bpstart bpstop bpisup bpparam register bpnworkers
#' @importFrom nipals nipals
#' @importFrom utils head
runNipalsSVD <- function(x, k=min(dim(x)), nu=k, nv=k, center=NULL, scale=NULL, deferred=FALSE, ..., fold=5L, BPPARAM=NULL)
{
    checked <- check_numbers(x, k=k, nu=nu, nv=nv)
    k <- checked$k
    nv <- checked$nv
    nu <- checked$nu

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
    FUN <- function(x, nu, nv, ...) { 
        # FUN is assumed to only take 'nu' and 'nv',
        # so we need to fill-in nipals()'s version of 'k'.
        out <- nipals(x, ncomp=max(nu, nv), ..., scale=FALSE, center=FALSE)
        colnames(out$scores) <- colnames(out$loadings) <- NULL
        list(u=out$scores, d=out$eig, v=out$loadings) 
    }

    if (use_crossprod(x, fold)) {
        res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=FUN, ...)
    } else {
        res <- FUN(x, nu=nu, nv=nv, ...)
        res$d <- head(res$d, k)
    }
    
    res
}

