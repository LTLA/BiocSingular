#' @importFrom BiocGenerics t
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMeans2 rowSums2
standardize_matrix <- function(x, center=NULL, scale=NULL) 
# Centers and scales the matrix, if this is required.
{
    if (is.null(center) && is.null(scale)) {
        return(x)
    }
    
    dx <- DelayedArray(x)
    dx <- t(dx)

    if (!is.null(center)) {
        if (is.logical(center)) {
            center <- rowMeans2(dx)
        }
        dx <- dx - center
    }

    if (!is.null(scale)) {
        if (is.logical(scale)) {
            scale <- sqrt(rowSums(dx^2) / (ncol(dx) - 1)) # consistent with 'scale' for non-default 'center'.
        }
        dx  <- dx/scale
    }
    
    t(dx)
}
