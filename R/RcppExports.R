# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @useDynLib BiocSingular
#' @importFrom Rcpp sourceCpp
compute_scale <- function(mat, centering) {
    .Call('_BiocSingular_compute_scale', PACKAGE = 'BiocSingular', mat, centering)
}

