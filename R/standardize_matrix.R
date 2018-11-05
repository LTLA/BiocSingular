#' @importFrom DelayedArray DelayedArray
standardize_matrix <- function(x, center=NULL, scale=NULL, deferred=FALSE)
# Creates a deferred or delayed centered and scaled matrix.
# The two choices have different implications for speed and accuracy.
{
    if (deferred) {
        X <- bs_matrix(x, center=center, scale=scale)
    } else {
        X <- as.matrix(x) # DelayedArray(x) # uncomment once DA supports crossprod, etc.
        if (!is.null(center)) {
            X <- sweep(X, 2, center, "-", check.margin=FALSE)
        }
        if (!is.null(scale)) {
            X <- sweep(X, 2, scale, "/", check.margin=FALSE)
        }
    }
    return(X)
}

