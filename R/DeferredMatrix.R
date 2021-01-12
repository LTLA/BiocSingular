#' @title The DeferredMatrix class
#'
#' @description 
#' This has been deprecated in favor of the \linkS4class{ScaledMatrix} class 
#' from the \pkg{ScaledMatrix} package - use those constructors instead.
#'
#' @docType class
#' @name DeferredMatrix
#' @aliases
#' DeferredMatrixSeed
NULL

#' @export
#' @importFrom ScaledMatrix ScaledMatrixSeed
DeferredMatrixSeed <- function(x, center=NULL, scale=NULL) {
    .Deprecated("ScaledMatrixSeed")
    ScaledMatrixSeed(x, center=center, scale=scale)
}

#' @export
#' @importFrom ScaledMatrix ScaledMatrix
DeferredMatrix <- function(x, center=NULL, scale=NULL) {
    .Deprecated("ScaledMatrix")
    ScaledMatrix(x, center=center, scale=scale)
}
