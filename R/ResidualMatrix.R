# BiocSingular matrix class for a reconstructed low-rank matrix.
# Useful for representing the low-rank matrix without actually making it.

###################################
###################################
###################################
# Constructing the seed.

#' @export
#' @importFrom methods new is
#' @importFrom Matrix crossprod
ResidualMatrixSeed <- function(x, design=matrix(1, ncol(x), 1)) {
    if (missing(x)) {
        x <- matrix(0, 0, 0)
    } else if (is(x, "ResidualMatrixSeed")) {
        return(x)
    } 

    QR <- qr(design)
    Q <- qr.Q(QR)
    new("ResidualMatrixSeed", .matrix=x, Q=Q, Qty=crossprod(Q, x), transposed=FALSE)
}

#' @importFrom S4Vectors setValidity2
#' @importFrom methods is
setValidity2("ResidualMatrixSeed", function(object) {
    msg <- character(0)

    x <- get_matrix2(object)
    Q <- get_Q(object)
    if (nrow(x)!=nrow(Q)) {
        msg <- c(msg, "'nrow(x)' and 'nrow(Q)' are not the same")
    }
    if (!is.numeric(Q)) {
        msg <- c(msg, "'Q' should be a numeric matrix")
    }

    Qty <- get_Qty(object)
    if (ncol(x)!=ncol(Qty)) {
        msg <- c(msg, "'ncol(x)' and 'ncol(Qty)' are not the same")
    }
    if (ncol(Q)!=nrow(Qty)) {
        msg <- c(msg, "'ncol(Q)' and 'nrow(Qty)' are not the same")
    }
    if (!is.numeric(Q)) {
        msg <- c(msg, "'Qty' should be a numeric matrix")
    }

    if (length(is_transposed(object))!=1L) {
        msg <- c(msg, "'transposed' must be a logical scalar")
    } 

    if (length(msg)) {
        return(msg)
    } 
    return(TRUE)
})

#' @export
#' @importFrom methods show
setMethod("show", "ResidualMatrixSeed", function(object) {
    cat(sprintf("%i x %i ResidualMatrixSeed object", nrow(object), ncol(object)),
    sep="\n")
})

###################################
# Internal getters.

get_Q <- function(x) x@Q

get_Qty <- function(x) x@Qty

###################################
# DelayedArray support utilities. 

#' @export
setMethod("dim", "ResidualMatrixSeed", function(x) {
    d <- dim(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
setMethod("dimnames", "ResidualMatrixSeed", function(x) {
    d <- dimnames(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
#' @importFrom DelayedArray extract_array
#' @importFrom Matrix t crossprod
setMethod("extract_array", "ResidualMatrixSeed", function(x, index) {
    if (was_transposed <- is_transposed(x)) {
        x <- transpose_ResidualMatrixSeed(x)
        index <- rev(index)
    }
	x2 <- subset_ResidualMatrixSeed(x, index[[1]], index[[2]])
    resid <- get_matrix2(x2) - get_Q(x2) %*% get_Qty(x2)
    if (was_transposed) {
        resid <- t(resid)
    }
    resid
})

###################################
# Additional utilities for efficiency.

subset_ResidualMatrixSeed <- function(x, i, j) {
    mat <- get_matrix2(x)
    Q <- get_Q(x)
    Qty <- get_Qty(x)

    if (!is.null(i)) {
        mat <- mat[i,,drop=FALSE]
        Q <- Q[i,,drop=FALSE]
    }
    if (!is.null(j)) {
        mat <- mat[,j,drop=FALSE]
        Qty <- Qty[,j,drop=FALSE]
    }

    initialize(x, .matrix=mat, Q=Q, Qty=Qty)
}

transpose_ResidualMatrixSeed <- function(x) {
    initialize(x, transposed=!is_transposed(x))
}

rename_ResidualMatrixSeed <- function(x, value) {
    x2 <- get_matrix2(x)
    dimnames(x2) <- value
    initialize(x, .matrix=x2)
}

###################################
###################################
###################################
# Constructing the matrix.

#' @export
#' @importFrom DelayedArray DelayedArray
ResidualMatrix <- function(x, design) {
    DelayedArray(ResidualMatrixSeed(x, design))
}

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "ResidualMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="ResidualMatrix")
)

###################################
# Overridden utilities from DelayedArray, for efficiency.

#' @export
#' @importFrom DelayedArray DelayedArray seed
setReplaceMethod("dimnames", "ResidualMatrix", function(x, value) {
    DelayedArray(rename_ResidualMatrixSeed(seed(x), value))
})

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("t", "ResidualMatrix", function(x) {
    DelayedArray(transpose_ResidualMatrixSeed(seed(x)))
})

#' @export
#' @importFrom DelayedArray DelayedArray seed
setMethod("[", "ResidualMatrix", function(x, i, j, ..., drop=TRUE) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL

    rseed <- seed(x)
    if (was_transposed <- is_transposed(rseed)) {
        rseed <- transpose_ResidualMatrixSeed(rseed)
        tmp <- i
        i <- j
        j <- tmp
    }

    rseed <- subset_ResidualMatrixSeed(rseed, i, j)
    if (was_transposed) {
        rseed <- transposedResidualMatrixSeed(rseed)
    }
    DelayedArray(rseed)
})

###################################
# Basic matrix stats.

#' @export
#' @importFrom Matrix colSums rowSums drop
setMethod("colSums", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(rowSums(t(x)))
    }

    out <- rep(1, nrow(x)) %*% x
    out <- drop(out)
    names(out) <- colnames(x)
    out
})

#' @export
#' @importFrom Matrix colSums rowSums drop
setMethod("rowSums", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(colSums(t(x)))
    }

    out <- x %*% rep(1, ncol(x))
    out <- drop(out)
    names(out) <- rownames(x)
    out
})

#' @export
#' @importFrom Matrix colMeans colSums
setMethod("colMeans", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) colSums(x)/nrow(x))

#' @export
#' @importFrom Matrix rowMeans rowSums
setMethod("rowMeans", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) rowSums(x)/ncol(x))

###################################
# Matrix multiplication.

# Note that "%*%" methods should NOT call other matrix multiplication methods
# directly on their ResidualMatrix arguments. Rather, any ResidualMatrix should
# be broken down into the seed or the underlying matrix before further multiplication.
# This reduces the risk of infinite S4 recursion when 'y' is also an S4 matrix class. 
# 
# Specifically, the .*_ResidualMatrix functions take a seed object that IGNORES
# any non-FALSE setting of @transposed. It will then break down the seed into 
# its constituents and perform the multiplication, such that any further S4
# dispatch occurs on the lower components of the seed. 
#
# We also coerce each matrix product to a full matrix - which it usually is, anyway 
# - to avoid the unnecessary overhead of multiplying DelayedArray instances.

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ResidualMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(t(y), x_seed))
    } else {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod
.rightmult_ResidualMatrix <- function(x_seed, y) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that ncol(y) is very small.
    as.matrix(get_matrix2(x_seed) %*% y) - get_Q(x_seed) %*% as.matrix(get_Qty(x_seed) %*% y)
}

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ANY", "ResidualMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        if (!is.null(dim(x))) {
            # Vectors don't quite behave as 1-column matrices here.
            # so we need to be a bit more careful.
            x <- t(x) 
        }
        out <- t(.rightmult_ResidualMatrix(y_seed, x))
    } else {
        out <- .leftmult_ResidualMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.leftmult_ResidualMatrix <- function(x, y_seed) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that nrow(x) is very small.
    as.matrix(x %*% get_matrix2(y_seed)) - as.matrix(x %*% get_Q(y_seed)) %*% get_Qty(y_seed)
}

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ResidualMatrix", "ResidualMatrix"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(t(y), x_seed))
    } else {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

###################################
# Crossproduct.

# Technically, this could be defined in terms of '%*%'.
# However, we re-define it in terms of 'crossprod' for efficiency,
# to exploit the potential availability of optimized crossprod for .matrix.

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ResidualMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        #  No need for t(), it's symmetric.
        out <- .tcp_ResidualMatrix(x_seed)
    } else {
        out <- .crossprod_ResidualMatrix(x_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod 
.crossprod_ResidualMatrix <- function(x_seed) {
    mat <- get_matrix2(x_seed)
    Qty <- get_Qty(x_seed)
    Q <- get_Q(x_seed)

    # We assume that nrow(Q) >> ncol(Q) and nrow(mat) >> ncol(mat) 
    # in order for this to be efficient.
    ytQQty <- crossprod(Qty)
    QtQ <- crossprod(Q)

    # Using this addition order to minimize numeric instability.
    (crossprod(mat) - ytQQty) + (crossprod(Qty, QtQ %*% Qty) - ytQQty)
}

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ResidualMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    } else {
        out <- .rightcross_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod
.rightcross_ResidualMatrix <- function(x_seed, y) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that ncol(y) is very small.
    as.matrix(crossprod(get_matrix2(x_seed), y)) - crossprod(get_Qty(x_seed), as.matrix(crossprod(get_Q(x_seed), y)))
}

#' @export
#' @importFrom Matrix crossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ANY", "ResidualMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        out <- t(.rightmult_ResidualMatrix(y_seed, x))
    } else {
        out <- .leftcross_ResidualMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod
.leftcross_ResidualMatrix <- function(x, y_seed) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that ncol(x) is very small.
    as.matrix(crossprod(x, get_matrix2(y_seed))) - as.matrix(crossprod(x, get_Q(y_seed))) %*% get_Qty(y_seed)
}

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ResidualMatrix", "ResidualMatrix"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    } else {
        out <- .rightcross_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

###################################
# Transposed crossproduct.

# Technically, this could be defined in terms of '%*%'.
# However, we re-define it in terms of 'tcrossprod' for efficiency,
# to exploit the potential availability of optimized crossprod for .matrix.

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ResidualMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .crossprod_ResidualMatrix(x_seed)
    } else {
        out <- .tcp_ResidualMatrix(x_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod 
.tcp_ResidualMatrix <- function(x_seed) {
    mat <- get_matrix2(x_seed)
    Qty <- get_Qty(x_seed)
    Q <- get_Q(x_seed)

    # We assume that ncol(mat) >> nrow(mat) for this to be efficient.
    # We also try to avoid constructing QQt under the assumption that 
    # nrow(Q) >> ncol(Q) for Q derived from a full-rank design matrix.
    QQtyyt<- Q %*% tcrossprod(Qty, mat)
    QQtyytQQt <- tcrossprod(QQtyyt %*% Q, Q)

    # Using this addition order to minimize numeric instability.
    (tcrossprod(mat) - QQtyyt) + (QQtyytQQt - t(QQtyyt))
}

#' @export
#' @importFrom Matrix tcrossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ResidualMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(y, x_seed))
    } else {
        out <- .righttcp_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.righttcp_ResidualMatrix <- function(x_seed, y) {
    # Order of operations chosen to minimize size of intermediates,
    # assuming that nrow(y) is very small.
    as.matrix(tcrossprod(get_matrix2(x_seed), y)) - get_Q(x_seed) %*% as.matrix(tcrossprod(get_Qty(x_seed), y))
}

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ANY", "ResidualMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        out <- .leftmult_ResidualMatrix(x, y_seed)
    } else {
        out <- .lefttcp_ResidualMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.lefttcp_ResidualMatrix <- function(x, y_seed) {
    # Order of operations chosen to minimize size of intermediates.
    # assuming that nrow(x) is very small.
    as.matrix(tcrossprod(x, get_matrix2(y_seed))) - tcrossprod(as.matrix(tcrossprod(x, get_Qty(y_seed))), get_Q(y_seed))
}

#' @export
#' @importFrom Matrix tcrossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ResidualMatrix", "ResidualMatrix"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(y, x_seed))
    } else {
        out <- .righttcp_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})
