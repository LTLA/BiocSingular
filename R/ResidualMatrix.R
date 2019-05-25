# BiocSingular matrix class for a reconstructed low-rank matrix.
# Useful for representing the low-rank matrix without actually making it.

###################################
###################################
###################################
# Constructing the seed.

#' @export
#' @importFrom methods new is
ResidualMatrixSeed <- function(x, design=matrix(1, ncol(x), 1)) {
    if (missing(x)) {
        x <- matrix(0, 0, 0)
    } else if (is(x, "ResidualMatrixSeed")) {
        return(x)
    } 

    QR <- qr(design)
    Q <- qr.Q(QR)
    new("ResidualMatrixSeed", .matrix=x, QR=QR, Q=Q, transposed=FALSE)
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

    QR <- get_QR(object)
    if (!is(QR, "qr")) {
        msg <- c(msg, "'QR' should be a list of class 'qr'")
    } 
    if (nrow(x)!=nrow(QR$qr)) {
        msg <- c(msg, "'nrow(QR$qr)' and 'nrow(x)' are not the same")
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

get_QR <- function(x) x@QR

get_Q <- function(x) x@Q

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
#' @importFrom Matrix t 
setMethod("extract_array", "ResidualMatrixSeed", function(x, index) {
    if (was_transposed <- is_transposed(x)) {
        x <- transpose_ResidualMatrixSeed(x)
        index <- rev(index)
    }
	x2 <- subset_ResidualMatrixSeed(x, index[[2]])

    resid <- qr.resid(get_QR(x2), get_matrix2(x2))
    dimnames(resid) <- dimnames(x2)
    if (!is.null(index[[1]])) {
        resid <- resid[index[[1]],,drop=FALSE]
    }
    if (was_transposed) {
        resid <- t(resid)
    }
    resid
})

###################################
# Additional utilities for efficiency.

subset_ResidualMatrixSeed <- function(x, j) {
    if (is.null(j)) {
        x
    } else {
        initialize(x, .matrix=get_matrix2(x)[,j,drop=FALSE])
    }
}

transpose_ResidualMatrixSeed <- function(x) {
    initialize(x, transposed=!get_transposed(x))
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

    rseed <- subset_ResidualMatrixSeed(rseed, j)
    x <- DelayedArray(rseed)

    if (is.null(i)) {
        if (was_transposed) {
            rseed <- tranpose_ResidualMatrixSeed(rseed)
        }
        x <- DelayedArray(rseed)
    } else {
        x <- DelayedArray(rseed)
        x <- callNextMethod(x=x, i=i) # Becomes a DelayedMatrix.
        if (was_transposed) {
            x <- t(x)
        }
    }
    x
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

#' @export
#' @importFrom Matrix crossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ResidualMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        return(t(t(y) %*% t(x)))    
    }

    X <- get_matrix2(x_seed)
    component1 <- X %*% y
    Q <- get_Q(x_seed)
    component2 <- Q %*% crossprod(Q, component1)

    DelayedArray(component1 - component2)
})

#' @export
#' @importFrom Matrix tcrossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ANY", "ResidualMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        return(t(t(y) %*% t(x)))    
    }

    Y <- get_matrix2(y_seed)
    component1 <- x %*% Y
    Q <- get_Q(y_seed)
    component2 <- tcrossprod(x %*% Q, Q) %*% Y

    DelayedArray(component1 - component2)
})
