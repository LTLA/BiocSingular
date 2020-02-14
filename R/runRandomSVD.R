#' @export
#' @importFrom BiocParallel SerialParam bpstart bpstop 
#' @importFrom utils head
#' @importClassesFrom DelayedArray DelayedMatrix
#' @importFrom methods is
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
runRandomSVD <- function(x, k=5, nu=k, nv=k, center=FALSE, scale=FALSE, deferred=FALSE, ..., fold=Inf, BPPARAM=SerialParam())
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
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    if (!.bpisup2(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred, BPPARAM=BPPARAM)

    if (use_crossprod(x, fold)) {
        FUN <- function(x, nu, nv, ...) { 
            # FUN is assumed to only take 'nu' and 'nv',
            # so we need to fill-in rsvd()'s version of 'k'.
            safe_rsvd(x, k=max(nu, nv), nu=nu, nv=nv, ...)
        }
        res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=FUN, ...) 

    } else {
        res <- safe_rsvd(x, k=max(nu, nv, k), nu=nu, nv=nv, ...) 
        res$d <- head(res$d, k)
        res <- standardize_output_SVD(res, x)
    }

    res
}

#' @importFrom rsvd rsvd 
safe_rsvd <- function(x, k, nu, nv, ...)
# Wrapper that guarantees return of U, D and V, even if nu or nv are zero.
{
    out <- rsvd(x, k=k, nu=nu, nv=nv, ...)
    if (!nu) {
        out$u <- matrix(0, nrow(x), 0)
    }
    if (!nv) {
        out$v <- matrix(0, ncol(x), 0)
    }
    out
}
