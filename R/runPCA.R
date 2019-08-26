#' Principal components analysis
#'
#' Perform a principal components analysis (PCA) on a target matrix with a specified SVD algorithm.
#' 
#' @param x A numeric matrix-like object with samples as rows and variables as columns.
#' @param rank Integer scalar specifying the number of principal components to retain.
#' @param center A logical scalar indicating whether columns of \code{x} should be centered before the PCA is performed.
#' Alternatively, a numeric vector of length \code{ncol(x)} containing the value to subtract from each column of \code{x}.
#' @param scale A logical scalar indicating whether columns of \code{x} should be scaled to unit variance before the PCA is performed.
#' Alternatively, a numeric vector of length \code{ncol(x)} containing the scaling factor for each column of \code{x}.
#' @param get.rotation A logical scalar indicating whether rotation vectors should be returned.
#' @param get.pcs A logical scalar indicating whether the principal component scores should be returned.
#' @param as.lem Logical scalar specifying whether the results should be returned as a \linkS4class{LinearEmbeddingMatrix}.
#' If \code{TRUE}, any \code{FALSE} values of \code{get.pcs} or \code{get.rotation} are ignored.
#' @param ... For the generic, additional arguments to pass to specific methods upon dispatch.
#' 
#' For the \code{ANY} method, further arguments to pass to \code{\link{runSVD}}.
#' This includes \code{BSPARAM} to specify the algorithm that should be used, and \code{BPPARAM} to control parallelization.
#' 
#' @details
#' This function simply calls \code{\link{runSVD}} and converts the results into a format similar to that returned by \code{\link{prcomp}}.
#'
#' The generic is exported to allow other packages to implement their own \code{runPCA} methods for other \code{x} objects, e.g., \pkg{scater} for SingleCellExperiment inputs.
#' 
#' @return
#' By default, a list is returned containing:
#' \itemize{
#' \item \code{sdev}, a numeric vector of length \code{rank} containing the standard deviations of the first \code{rank} principal components.
#' \item \code{rotation}, a numeric matrix with \code{rank} columns and \code{nrow(x)} rows, containing the first \code{rank} rotation vectors.
#' This is only returned if \code{get.rotation=TRUE}.
#' \item \code{x}, a numeric matrix with \code{rank} columns and \code{ncol(x)} rows, containing the scores for the first \code{rank} principal components.
#' This is only returned if \code{get.pcs=TRUE}.
#' }
#'
#' If \code{as.lem=TRUE}, a \linkS4class{LinearEmbeddingMatrix} is returned
#' containing the contents of \code{rotation} as the \code{\link{featureLoadings}};
#' the contents of \code{x} as the \code{\link{sampleFactors}};
#' and a \code{sdev} field in the \code{\link{factorData}}.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{runSVD}} for the underlying SVD function.
#' 
#' \code{?\linkS4class{BiocSingularParam}} for details on the algorithm choices.
#'
#' @examples
#' A <- matrix(rnorm(100000), ncol=20)
#' str(out <- runPCA(A, rank=10))
#'
#' # Or, as a LinearEmbeddingMatrix.
#' (out2 <- runPCA(A, rank=10, as.lem=TRUE))
#'
#' @name runPCA
NULL

#' @export
#' @rdname runPCA
#' @importFrom S4Vectors DataFrame
setMethod("runPCA", "ANY", function(x, rank, center=TRUE, scale=FALSE, get.rotation=TRUE, get.pcs=TRUE, as.lem=FALSE, ...) 
# Converts SVD results to PCA results.
{
    if (as.lem) {
        get.rotation <- get.pcs <- TRUE
    }

    svd.out <- runSVD(x, k=rank, 
        nu=ifelse(get.pcs, rank, 0),
        nv=ifelse(get.rotation, rank, 0),
        center=center, scale=scale, ...)

    out <- list(sdev=svd.out$d / sqrt(nrow(x) - 1))
    if (get.rotation) {
        out$rotation <- svd.out$v
        colnames(out$rotation) <- sprintf("PC%i", seq_len(ncol(out$rotation)))
    } 
    if (get.pcs) {
        out$x <- sweep(svd.out$u, 2, svd.out$d, "*")
        colnames(out$x) <- sprintf("PC%i", seq_len(ncol(out$x)))
    }

    if (as.lem) {
        out <- LinearEmbeddingMatrix(sampleFactors=out$x, 
            featureLoadings=out$rotation, 
            factorData=DataFrame(sdev=out$sdev))
    }

    out
})

#' @export
setMethod("runPCA", "ResidualMatrix", function(x, center=TRUE, ...) {
    if (center && is_centered(seed(x))) {
        center <- FALSE
    }
    callNextMethod()
})
