#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom irlba irlba
runIrlba <- function(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, extra.work=7, maxit=1000, tol=1e-5, BPPARAM=SerialParam())
# Wrapper for irlba(), switching to the appropriate multiplication algorithm for  
{
    if (nu==0 && nv==0 && k==0) {
        return(list(d=numeric(0),
                    u=matrix(0, nrow(x), 0),
                    v=matrix(0, ncol(x), 0)))
    }

    args <- list(A=x, nu=nu, nv=max(nv, k),
            work=max(k, nu, nv) + extra.work, 
            center=center, scale=scale, 
            maxit=maxit, tol=tol)

    # Fix when DA parallel multiplication is operational.
    if (bpnworkers(BPPARAM)!=1L) {
        args$mult <- `%*%`
    }

    res <- do.call(irlba, args)
    res$v <- res$v[,seq_len(nv),drop=FALSE]
    res[c("d", "u", "v")]
}
