# BiocSingular matrix class for a reconstructed low-rank matrix.
# Useful for representing the low-rank matrix without actually making it.

###################################
###################################
###################################
# Constructing the seed.

#' @export
ResidualMatrixSeed <- function(x, design=NULL) {
    .Deprecated(old="BiocSingular::ResidualMatrixSeed", new="ResidualMatrix::ResidualMatrixSeed")
    ResidualMatrix::ResidualMatrixSeed(x, design=design)
}

#' @export
#' @importFrom DelayedArray DelayedArray
ResidualMatrix <- function(x, design=NULL) {
    .Deprecated(old="BiocSingular::ResidualMatrix", new="ResidualMatrix::ResidualMatrix")
    DelayedArray(ResidualMatrix::ResidualMatrixSeed(x, design))
}

