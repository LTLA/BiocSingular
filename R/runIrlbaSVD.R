#' @export
#' @importFrom BiocParallel bpstart bpstop bpisup bpparam register
#' @importFrom irlba irlba
#' @importFrom utils head
#' @importClassesFrom Matrix dgCMatrix
runIrlbaSVD <- function(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, extra.work=7, ..., fold=5L, BPPARAM=NULL)
# Wrapper for irlba(), switching to the appropriate multiplication algorithm for  
{
    if (nu==0 && nv==0 && k==0) {
        return(list(d=numeric(0),
                    u=matrix(0, nrow(x), 0),
                    v=matrix(0, ncol(x), 0)))
    }

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

    args <- list(work=max(k, nu, nv) + extra.work, ...)

    if (use_crossprod(x, fold)) {
        x <- bs_matrix(x, center=center, scale=scale)
        res <- do.call(svd_via_crossprod, c(list(x, k=k, nu=nu, nv=nv, FUN=irlba, BPPARAM=BPPARAM), args))

    } else {
        args$nu <- nu
        args$nv <- max(k, nv)

        if (bpnworkers(BPPARAM)==1L && (is.matrix(x) || is(x, "dgCMatrix"))) {
            # Try to use 'fastpath=TRUE' if possible.
            res <- do.call(irlba, c(list(A=x, center=center, scale=scale), args))

        } else {
            x <- bs_matrix(x, center=center, scale=scale)
            res <- do.call(irlba, c(list(A=x, fastpath=FALSE), args))
        }

        res$v <- res$v[,seq_len(nv),drop=FALSE]
        res$d <- head(res$d, k)
    }

    res[c("d", "u", "v")]
}
