check_numbers <- function(x, k, nv, nu) 
# Sanity checks for 'k', 'nv' and 'nu'.
{
    if (!is.finite(k) || !is.finite(nv) || !is.finite(nu) || k < 0 || nv < 0 || nu < 0) {
        stop("number of requested singular values/vectors must be non-negative")
    }

    limit <- min(dim(x))
    if (max(k, nu, nv) > limit) {
        warning("more singular values/vectors requested than available")
        k <- min(limit, k)
        nv <- min(limit, nv)
        nu <- min(limit, nu)
    }
    
    list(k=k, nv=nv, nu=nu)
}
