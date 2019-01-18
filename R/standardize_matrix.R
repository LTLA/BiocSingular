#' @importFrom DelayedArray DelayedArray sweep
standardize_matrix <- function(x, center=NULL, scale=NULL, deferred=FALSE)
# Creates a deferred or delayed centered and scaled matrix.
# The two choices have different implications for speed and accuracy.
{
    if (deferred) {
        X <- bs_matrix(x, center=center, scale=scale)
    } else {
        X <- DelayedArray(x)
        if (!is.null(center)) {
            X <- sweep(X, 2, center, "-") 
        }
        if (!is.null(scale)) {
            X <- sweep(X, 2, scale, "/")
        }
    }
    return(X)
}

standardize_output_SVD <- function(res) 
# Provide a common standard output for all SVD functions.
{
    res$d <- as.numeric(res$d)
    res$u <- as.matrix(res$u)
    res$v <- as.matrix(res$v)
    res[c("d", "u", "v")]
}
