#' @export
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup
#' @importFrom rsvd rsvd 
#' @importFrom utils head
#' @importClassesFrom DelayedArray DelayedMatrix
#' @importFrom methods is
runRandomSVD <- function(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, deferred=FALSE, ..., fold=5L, BPPARAM=NULL)
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

    x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred)

    if (use_crossprod(x, fold)) {
        FUN <- function(x, nu, nv, ...) { 
            # FUN is assumed to only take 'nu' and 'nv',
            # so we need to fill-in rsvd()'s version of 'k'.
            rsvd(x, k=max(nu, nv), nu=nu, nv=nv, ...)
        }
        res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=FUN, ...) 

    } else {
        res <- rsvd(x, k=max(nu, nv, k), nu=nu, nv=nv, ...) 
        res$d <- head(res$d, k)
    }

    res[c("d", "u", "v")]
}

