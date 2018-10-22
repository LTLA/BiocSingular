#' @export
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup
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

    args <- list(work=max(k, nu, nv) + extra.work, ...)

    if (use_crossprod(x, fold)) {
        x <- standardize_matrix(x, center=center, scale=scale)
        x <- as.matrix(x) # remove once crossprod supports DAs.
        res <- do.call(svd_via_crossprod, c(list(x, k=k, nu=nu, nv=nv, FUN=irlba, BPPARAM=BPPARAM), args))

    } else {
        if (bpnworkers(BPPARAM)!=1L) {
            if (!bpisup(BPPARAM)) {
                bpstart(BPPARAM) # lots of multiplications, so we set up the backend.
                on.exit(bpstop(BPPARAM))
            }
            x <- bpmatrix(x, BPPARAM)
        }

        res <- do.call(irlba, c(list(A=x, nu=nu, nv=max(k, nv), center=center, scale=scale), args))
        res$v <- res$v[,seq_len(nv),drop=FALSE]
        res$d <- head(res$d, k)
    }

    res[c("d", "u", "v")]
}
