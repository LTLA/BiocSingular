#' @export
#' @importFrom BiocParallel bpstart bpstop bpisup bpparam register bpnworkers SerialParam
#' @importFrom irlba irlba
#' @importFrom utils head
#' @importClassesFrom Matrix dgCMatrix
runIrlbaSVD <- function(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, deferred=FALSE, extra.work=7, ..., fold=5, BPPARAM=SerialParam())
# Wrapper for irlba(), switching to the appropriate multiplication algorithm for  
{
    if (nu==0 && nv==0 && k==0) {
        return(list(d=numeric(0),
                    u=matrix(0, nrow(x), 0),
                    v=matrix(0, ncol(x), 0)))
    }

    checked <- check_numbers(x, k=k, nu=nu, nv=nv)
    k <- checked$k
    nv <- checked$nv
    nu <- checked$nu

    # Setting up the parallelization environment.
    old <- bpparam()
    register(BPPARAM)
    on.exit(register(old))

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    args <- list(work=max(k, nu, nv) + extra.work, ...)

    if (use_crossprod(x, fold)) {
        x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred)
        res <- do.call(svd_via_crossprod, c(list(x, k=k, nu=nu, nv=nv, FUN=irlba, BPPARAM=BPPARAM), args))

    } else {
        args$nu <- nu
        args$nv <- max(k, nv)

        if (bpnworkers(BPPARAM)==1L && (is.matrix(x) || is(x, "dgCMatrix"))) {
            # Try to use 'fastpath=TRUE' if possible.
            res <- do.call(irlba, c(list(A=x, center=center, scale=scale), args))

        } else {
            x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred)
            res <- do.call(irlba, c(list(A=x, fastpath=FALSE), args))
        }

        res$v <- res$v[,seq_len(nv),drop=FALSE]
        res$d <- head(res$d, k)
        res <- standardize_output_SVD(res)
    }

    res
}
