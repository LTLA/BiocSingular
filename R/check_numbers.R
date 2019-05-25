check_numbers <- function(x, k, nv, nu, limit=min(dim(x))) 
# Sanity checks for 'k', 'nv' and 'nu'.
{
    if (k < 0 || nv < 0 || nu < 0) {
        stop("number of requested singular values/vectors must be non-negative")
    }

    if (is.infinite(k)) { 
        k <- limit
    }
    if (is.infinite(nv)) {
        nv <- limit
    }
    if (is.infinite(nu)) {
        nu <- limit
    }

    if (max(k, nu, nv) > limit) {
        warning("more singular values/vectors requested than available")
        k <- min(limit, k)
        nv <- min(limit, nv)
        nu <- min(limit, nu)
    }
    
    list(k=k, nv=nv, nu=nu)
}
