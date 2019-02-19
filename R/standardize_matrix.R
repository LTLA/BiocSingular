#' @importFrom DelayedArray DelayedArray 
#' @importFrom BiocParallel SerialParam bpnworkers
standardize_matrix <- function(x, center=NULL, scale=NULL, deferred=FALSE, BPPARAM=SerialParam())
# Creates a deferred or delayed centered and scaled matrix.
# The two choices have different implications for speed and accuracy.
{
    stats <- .compute_center_and_scale(x, center, scale)
    center <- stats$center
    scale <- stats$scale

    if (deferred) {
        if (bpnworkers(BPPARAM)==1L) {
            original <- x # exploit original (non-DA) matrix mult functions.
        } else {
            original <- DelayedArray(x) # exploit parallelization for DAs.
        }
        X <- DeferredMatrix(original, center=center, scale=scale)

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

#' @importFrom DelayedArray DelayedArray sweep
#' @importFrom BiocGenerics colSums 
.compute_center_and_scale <- function(x, center, scale) {
    if (is.logical(center)) {
        if (center) {
            center <- .safe_colSums(x)/nrow(x)
        } else {
            center <- NULL
        }
    }

    if (is.logical(scale)) {
        if (scale) {
            sub <- DelayedArray(x)
            if (!is.null(center)) {
                sub <- sweep(x, 2, center, "-")
            }
            scale <- colSums(sub^2) / (nrow(x) - 1L) # mimic scale() behaviour for any 'center'.
            scale <- sqrt(scale)
        } else {
            scale <- NULL
        }
    }

    list(center=center, scale=scale) 
}

standardize_output_SVD <- function(res) 
# Provide a common standard output for all SVD functions.
{
    res$d <- as.numeric(res$d)
    res$u <- as.matrix(res$u)
    res$v <- as.matrix(res$v)
    res[c("d", "u", "v")]
}
