#' @importFrom DelayedArray DelayedArray 
#' @importFrom ScaledMatrix ScaledMatrix
#' @importFrom BiocParallel SerialParam bpnworkers
standardize_matrix <- function(x, center=FALSE, scale=FALSE, deferred=FALSE, BPPARAM=SerialParam())
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
        X <- ScaledMatrix(original, center=center, scale=scale)

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

#' @importFrom Matrix colMeans
#' @importFrom beachmat colBlockApply
.compute_center_and_scale <- function(x, center, scale) {
    if (is.logical(center)) {
        if (center) {
            center <- colMeans(x)
        } else {
            center <- NULL
        }
    }

    if (is.logical(scale)) {
        if (scale) {
            # mimic scale() behaviour for any 'center'.
            scale <- colBlockApply(x, FUN=.compute_scale, center=center)
            scale <- unlist(scale)
        } else {
            scale <- NULL
        }
    }

    list(center=center, scale=scale) 
}

#' @importFrom DelayedArray makeNindexFromArrayViewport currentViewport
.compute_scale <- function(block, center) {
    if (!is.null(center)) {
        vp <- currentViewport()
        cols <- makeNindexFromArrayViewport(vp, expand.RangeNSBS=TRUE)[[2]]
        if (!is.null(cols)) {
            center <- center[cols]
        }
    }
    compute_scale(block, center)
}

standardize_output_SVD <- function(res, x) 
# Provide a common standard output for all SVD functions.
{
    res$d <- as.numeric(res$d)
    res$u <- as.matrix(res$u)
    rownames(res$u) <- rownames(x)
    res$v <- as.matrix(res$v)
    rownames(res$v) <- colnames(x)
    res[c("d", "u", "v")]
}
