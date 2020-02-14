#' @export
#' @importFrom BiocParallel bpstart bpstop bpnworkers SerialParam
#' @importFrom irlba irlba
#' @importFrom utils head
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
runIrlbaSVD <- function(x, k=5, nu=k, nv=k, center=FALSE, scale=FALSE, deferred=FALSE, extra.work=7, ..., 
    fold=Inf, BPPARAM=SerialParam())
{
    if (nu==0 && nv==0 && k==0) {
        return(list(d=numeric(0),
                    u=matrix(0, nrow(x), 0),
                    v=matrix(0, ncol(x), 0)))
    }

    # irlba() errors when asking for exactly the maximum number of SVs,
    # hence the lowered 'limit'.
    checked <- check_numbers(k=k, nu=nu, nv=nv, limit=min(dim(x))-1L)
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

    args <- list(work=max(k, nu, nv) + extra.work, ...)

    if (use_crossprod(x, fold)) {
        x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred, BPPARAM=BPPARAM)
        res <- do.call(svd_via_crossprod, c(list(x, k=k, nu=nu, nv=nv, FUN=irlba, BPPARAM=BPPARAM), args))

    } else {
        args$nu <- nu
        args$nv <- max(k, nv)

        if (bpnworkers(BPPARAM)==1L) {
            stats <- .compute_center_and_scale(x, center, scale)

            # Use irlba's native 'center' and 'scale', avoid overhead of Delayed/DeferredMatrix wrapper in 'standardize_matrix'.
            # Also try to use 'fastpath=TRUE' if possible, depending on the class of 'x'.
            res <- do.call(irlba, c(args, list(A=x, center=stats$center, scale=stats$scale)))

        } else {
            x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred, BPPARAM=BPPARAM)
            res <- do.call(irlba, c(list(A=x, fastpath=FALSE), args))
        }

        res$v <- res$v[,seq_len(nv),drop=FALSE]
        res$d <- head(res$d, k)
        res <- standardize_output_SVD(res, x)
    }

    res
}
