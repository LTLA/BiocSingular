# Wrapper class that enables parallelized matrix multiplication
# in IRLBA() and RSVD(), without needing to define a custom 
# matrix multiplication function.

#' @import methods
#' @importClassesFrom BiocParallel BiocParallelParam
setClass("bpmatrix", slots=c(MATRIX="ANY", BPPARAM="BiocParallelParam"))

get_matrix <- function(x) x@MATRIX

get_bpparam <- function(x) x@BPPARAM

#' @importFrom methods new
bpmatrix <- function(x, BPPARAM) new("bpmatrix", MATRIX=x, BPPARAM=BPPARAM) 

# Utilities.
setMethod("dim", "bpmatrix", function(x) dim(get_matrix(x)))

setMethod("dimnames", "bpmatrix", function(x) dimnames(get_matrix(x)))

setMethod("length", "bpmatrix", function(x) length(get_matrix(x)))

# Matrix multiplication.
setMethod("%*%", c("bpmatrix", "ANY"),
    function(x, y) bpmult(get_matrix(x), y, BPPARAM=get_bpparam(x))
)

setMethod("%*%", c("ANY", "bpmatrix"),
    function(x, y) bpmult(x, get_matrix(y), BPPARAM=get_bpparam(y))
)

setMethod("%*%", c("bpmatrix", "bpmatrix"),
    function(x, y) bpmult(get_matrix(x), get_matrix(y), BPPARAM=get_bpparam(x))
)

# Cross-product.
setMethod("crossprod", c("bpmatrix", "missing"),
    function(x, ...) bpcross(get_matrix(x), BPPARAM=get_bpparam(x))
)

setMethod("crossprod", c("bpmatrix", "ANY"),
    function(x, y, ...) bpcross(get_matrix(x), y, BPPARAM=get_bpparam(x))
)

setMethod("crossprod", c("ANY", "bpmatrix"),
    function(x, y, ...) bpcross(x, get_matrix(y), BPPARAM=get_bpparam(y))
)

setMethod("crossprod", c("bpmatrix", "bpmatrix"),
    function(x, y, ...) bpcross(get_matrix(x), get_matrix(y), BPPARAM=get_bpparam(x))
)

# Transposed cross-product.
setMethod("tcrossprod", c("bpmatrix", "missing"),
    function(x, ...) bptcross(get_matrix(x), BPPARAM=get_bpparam(x))
)

setMethod("tcrossprod", c("bpmatrix", "ANY"),
    function(x, y, ...) bptcross(get_matrix(x), y, BPPARAM=get_bpparam(x))
)

setMethod("tcrossprod", c("ANY", "bpmatrix"),
    function(x, y, ...) bptcross(x, get_matrix(y), BPPARAM=get_bpparam(y))
)

setMethod("tcrossprod", c("bpmatrix", "bpmatrix"),
    function(x, y, ...) bptcross(get_matrix(x), get_matrix(y), BPPARAM=get_bpparam(x))
)
