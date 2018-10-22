#' @export
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup
#' @importFrom rsvd rsvd 
#' @importFrom utils head
runRandomSVD <- function(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, extra.work=7, ..., fold=5L, BPPARAM=SerialParam())
# Wrapper for irlba(), switching to the appropriate multiplication algorithm for  
{
    if (nu==0 && nv==0 && k==0) {
        return(list(d=numeric(0),
                    u=matrix(0, nrow(x), 0),
                    v=matrix(0, ncol(x), 0)))
    }
        
    x <- standardize_matrix(x, center=center, scale=scale)

    if (use_crossprod(x, fold)) {
        x <- as.matrix(x) # remove once crossprod supports DAs.

        FUN <- function(x, nu, nv, ...) { 
            # FUN is assumed to only take 'nu' and 'nv',
            # so we need to fill-in rsvd()'s version of 'k'.
            rsvd(x, k=max(nu, nv), nu=nu, nv=nv, ...)
        }

        res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=FUN, ..., BPPARAM=BPPARAM)

    } else {
        if (bpnworkers(BPPARAM)!=1L) {
            if (!bpisup(BPPARAM)) {
                bpstart(BPPARAM) # lots of multiplications, so we set up the backend.
                on.exit(bpstop(BPPARAM))
            }
            x <- bpmatrix(x, BPPARAM)
        }

        res <- rsvd(x, k=max(nu, nv, k), nu=nu, nv=nv, ...) 
        res$d <- head(res$d, k)
    }

    res[c("d", "u", "v")]
}

