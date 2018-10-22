#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom irlba irlba
#' @importFrom utils head
runIrlba <- function(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, extra.work=7, ..., fold=5L, BPPARAM=SerialParam())
# Wrapper for irlba(), switching to the appropriate multiplication algorithm for  
{
    if (nu==0 && nv==0 && k==0) {
        return(list(d=numeric(0),
                    u=matrix(0, nrow(x), 0),
                    v=matrix(0, ncol(x), 0)))
    }

    args <- list(A=x, nu=nu, nv=max(nv, k),
            work=max(k, nu, nv) + extra.work, ...)

    if (bpnworkers(BPPARAM)!=1L) {
        args$mult <- function(x, y) { as.matrix(bpmult(x, y, BPPARAM=BPPARAM)) }
        args$fastpath <- FALSE
    }

    if (use_crossprod(x, fold)) {
        x <- standardize_matrix(x, center=center, scale=scale)
        x <- as.matrix(x) # remove once crossprod supports DAs.
        res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=safe_svd, BPPARAM=BPPARAM)
    } else {
        args$center <- center
        args$scale <- scale
        res <- do.call(irlba, args)
        res$v <- res$v[,seq_len(nv),drop=FALSE]
        res$d <- head(res$d, k)
    }

    res[c("d", "u", "v")]
}
