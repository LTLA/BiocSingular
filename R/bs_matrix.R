# BiocSingular internal matrix with better support for centering and scaling.
# Specifically, we multiply first and apply the centering/scaling afterwards.

#' @import methods
setClass("bs_matrix", slots=c(.matrix="ANY", center="numeric", scale="numeric", use_center="logical", use_scale="logical"))

#' @importFrom methods new
bs_matrix <- function(x, center=NULL, scale=NULL) {
    use_center <- !is.null(center)
    if (use_center && length(center)!=ncol(x)) {
        stop("length of 'center' must equal 'ncol(x)'")
    }
    
    use_scale <- !is.null(scale)
    if (use_scale && length(scale)!=ncol(x)) {
        stop("length of 'scale' must equal 'ncol(x)'")
    }

    new("bs_matrix", .matrix=x, center=as.numeric(center), scale=as.numeric(scale), use_center=use_center, use_scale=use_scale)
}

###################################
# Getters. 

get_matrix2 <- function(x) x@.matrix

get_center <- function(x) x@center

get_scale <- function(x) x@scale

use_center <- function(x) x@use_center

use_scale <- function(x) x@use_scale

###################################
# Utilities.

setMethod("dim", "bs_matrix", function(x) dim(get_matrix2(x)))

setMethod("dimnames", "bs_matrix", function(x) dimnames(get_matrix2(x)))

setMethod("length", "bs_matrix", function(x) length(get_matrix2(x)))

###################################
# Matrix multiplication.

# We attempt to use the '%*%' defined for '.matrix' in the 'bs_matrix'.
# This avoids expensive modifications such as loss of sparsity.
# Centering and scaling are factored out into separate operations.
#
# We assume that the non-'bs_matrix' argument is small and can be modified cheaply.
# We also assume that the matrix product is small and can be modified cheaply.
# This allows centering and scaling to be applied *after* multiplication.

setMethod("%*%", c("bs_matrix", "ANY"), function(x, y) {
    if (use_scale(x)) {
        y <- y / get_scale(x) 
    }

    out <- as.matrix(get_matrix2(x) %*% y)

    if (use_center(x)) {
        out <- sweep(out, 2, as.numeric(get_center(x) %*% y), "-", check.margin=FALSE)
    }

    out
})

#' @importFrom BiocGenerics rowSums
setMethod("%*%", c("ANY", "bs_matrix"), function(x, y) {
    out <- as.matrix(x %*% get_matrix2(y))

    if (use_center(y)) {
        if (is.null(dim(x))) {
            out <- out - get_center(y) * sum(x)
        } else {
            out <- out - outer(rowSums(x), get_center(y), "*")
        }
    }

    if (use_scale(y)) {
        out <- sweep(out, 2, get_scale(y), "/", check.margin=FALSE)
    }

    out
})

setMethod("%*%", c("bs_matrix", "bs_matrix"), function(x, y) {
    stop("multiplication of two 'bs_matrix' objects is not yet supported")
})

###################################
# Cross-product. 

# Again, we attempt to use the 'crossprod' defined for '.matrix' in the 'bs_matrix'.

#' @importFrom BiocGenerics colSums
setMethod("crossprod", c("bs_matrix", "missing"), function(x, y) {
    x0 <- get_matrix2(x)
    out <- as.matrix(crossprod(x0))

    if (use_center(x)) {
        centering <- get_center(x)
        colsums <- colSums(x0)

        # Minus, then add, then minus, to mitigate cancellation.
        out <- out - outer(centering, colsums)
        out <- out + outer(centering, centering) * nrow(x0)
        out <- out - outer(colsums, centering)
    }

    if (use_scale(x)) {
        out <- sweep(out / get_scale(x), 2, get_scale(x), "/", check.margin=FALSE)
    }

    out
})

setMethod("crossprod", c("bs_matrix", "ANY"), function(x, y) {

})

setMethod("crossprod", c("ANY", "bs_matrix"), function(x, y) {

})

###################################
# Transposed cross-product. 

setMethod("tcrossprod", c("bs_matrix", "missing"), function(x, y) {
    new.x <- get_matrix2(x)
    if (use_scale(x)) {
        new.x <- sweep(new.x, 2, get_scale(x), "/", check.margin=FALSE) # Can't avoid this.
    }

    out <- as.matrix(tcrossprod(new.x))
    
    if (use_center(x)) {
        centering <- get_center(x)
        if (use_scale(x)) {
            centering <- centering / get_scale(x)
        }

        # Minus, then add, then minus, to mitigate cancellation.
        out <- sweep(out, 2, as.numeric(tcrossprod(centering, new.x)), "-", check.margin=FALSE)
        out <- out + sum(centering^2)
        out <- out - as.numeric(new.x %*% centering)
    }

    out
})
