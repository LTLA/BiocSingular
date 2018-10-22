#' @importFrom BiocGenerics t
#' @importFrom DelayedArray DelayedArray
standardize_matrix <- function(x, center=NULL, scale=NULL) 
# Centers and scales the matrix, if this is required.
{
    if (is.null(center) && is.null(scale)) {
        return(x)
    }
    dx <- DelayedArray(x)
    dx <- t(dx)

    if (!is.null(center)) {
        dx <- dx - center
    }
    if (!is.null(scale)) {
        dx  <- dx/scale
    }
    
    t(dx)
}
