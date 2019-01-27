# BiocSingular matrix class with support for deferred centering and scaling.
# Specifically, we multiply first and apply the centering/scaling afterwards.

#' @export
#' @import methods
setClass("DeferredMatrix", slots=c(.matrix="ANY", center="numeric", scale="numeric", use_center="logical", use_scale="logical", transposed="logical"))

#' @export
#' @importFrom methods new
DeferredMatrix <- function(x, center=NULL, scale=NULL) {
    use_center <- !is.null(center)
    use_scale <- !is.null(scale)
    new("DeferredMatrix", .matrix=x, center=as.numeric(center), scale=as.numeric(scale), use_center=use_center, use_scale=use_scale, transposed=FALSE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("DeferredMatrix", function(object) {
    msg <- character(0)

    # Checking scalars.
    if (length(use_center(object))!=1L) {
        msg <- c(msg, "'use_center' must be a logical scalar")
    } 
    if (length(use_scale(object))!=1L) {
        msg <- c(msg, "'use_scale' must be a logical scalar")
    } 
    if (length(is_transposed(object))!=1L) {
        msg <- c(msg, "'transposed' must be a logical scalar")
    } 

    # Checking vectors.
    if (use_center(object) && length(get_center(object))!=ncol(object)) {
        msg <- c(msg, "length of 'center' must equal 'ncol(object)'")
    }
    if (use_scale(object) && length(get_scale(object))!=ncol(object)) {
        msg <- c(msg, "length of 'scale' must equal 'ncol(object)'")
    }

    if (length(msg)) {
        return(msg)
    } 
    return(TRUE)
})

#' @export
#' @importFrom methods show
setMethod("show", "DeferredMatrix", function(object) {
    cat(sprintf("%i x %i DeferredMatrix object", nrow(object), ncol(object)),
        sprintf("representation: %s", class(get_matrix2(object))),
        sprintf("centering: %s", if (use_center(object)) "yes" else "no"),
        sprintf("scaling: %s", if (use_scale(object)) "yes" else "no"),
    sep="\n")
})

###################################
# Getters. 

get_matrix2 <- function(x) x@.matrix

get_center <- function(x) x@center

get_scale <- function(x) x@scale

use_center <- function(x) x@use_center

use_scale <- function(x) x@use_scale

is_transposed <- function(x) x@transposed

###################################
# Utilities.

#' @export
setMethod("dim", "DeferredMatrix", function(x) {
    d <- dim(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
setMethod("dimnames", "DeferredMatrix", function(x) {
    d <- dimnames(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
setMethod("length", "DeferredMatrix", function(x) length(get_matrix2(x)))

#' @export
#' @importFrom methods new
setMethod("t", "DeferredMatrix", function(x) {
    x@transposed <- !is_transposed(x)
    x
})

#' @export
#' @importFrom BiocGenerics t
setMethod("as.matrix", "DeferredMatrix", function(x) {
    out <- get_matrix2(x)
    if (use_scale(x) || use_center(x)) {
        out <- t(out)
        if (use_center(x)) {
            out <- out - get_center(x)
        }
        if (use_scale(x)) {
            out <- out / get_scale(x)
        }
        if (!is_transposed(x)) out <- t(out)
    } else {
        if (is_transposed(x)) out <- t(out) 
    }

    as.matrix(out)
})

###################################
# Matrix subset.

#' @export
setMethod("[", c(x="DeferredMatrix", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ..., drop=TRUE) {
    if (is_transposed(x)) {
        x <- t(x)
        if (!missing(i)) {
            x <- x[,i,drop=FALSE]
        }
        if (!missing(j)) {
            x <- x[j,,drop=FALSE]
        }
        x <- t(x)

    } else {
        if (!missing(i)) {
            x@.matrix <- get_matrix2(x)[i,,drop=FALSE]
        }
            
        if (!missing(j)) {
            x@.matrix <- get_matrix2(x)[,j,drop=FALSE]
            if (is.character(j)) {
                j <- match(j, colnames(x))
            }
            
            if (use_scale(x)) {
                x@scale <- get_scale(x)[j]
            }
    
            if (use_center(x)) {
                x@center <- get_center(x)[j]
            }
        }
    }

    if (drop && any(dim(x)==1L)) {
        return(drop(as.matrix(x)))
    }
    return(x)
})

###################################
# Matrix multiplication.

# We attempt to use operators defined for '.matrix' in the 'DeferredMatrix'.
# This avoids expensive modifications such as loss of sparsity.
# Centering and scaling are factored out into separate operations.
#
# We assume that the non-'DeferredMatrix' argument is small and can be modified cheaply.
# We also assume that the matrix product is small and can be modified cheaply.
# This allows centering and scaling to be applied *after* multiplication.

#' @export
#' @importFrom BiocGenerics t
setMethod("%*%", c("DeferredMatrix", "ANY"), function(x, y) {
    if (is_transposed(x)) {
        return(t(t(y) %*% t(x)))    
    }

    if (use_scale(x)) {
        y <- y / get_scale(x) 
    }

    out <- as.matrix(get_matrix2(x) %*% y)

    if (use_center(x)) {
        out <- sweep(out, 2, as.numeric(get_center(x) %*% y), "-", check.margin=FALSE)
    }

    out
})

#' @export
#' @importFrom BiocGenerics t
setMethod("%*%", c("ANY", "DeferredMatrix"), function(x, y) {
    if (is_transposed(y)) {
        if (!is.null(dim(x))) x <- t(x) # as vectors don't quite behave as 1-column matrices here.
        return(t(t(y) %*% x))    
    }

    out <- as.matrix(x %*% get_matrix2(y))

    if (use_center(y)) {
        if (is.null(dim(x))) {
            out <- out - get_center(y) * sum(x)
        } else {
            out <- out - outer(Matrix::rowSums(x), get_center(y), "*")
        }
    }

    if (use_scale(y)) {
        out <- sweep(out, 2, get_scale(y), "/", check.margin=FALSE)
    }

    out
})

setMethod("%*%", c("DeferredMatrix", "DeferredMatrix"), function(x, y) {
    stop("multiplication of two 'DeferredMatrix' objects is not yet supported")
})

###################################
# Cross-product. 

#' @export
#' @importFrom Matrix crossprod 
setMethod("crossprod", c("DeferredMatrix", "missing"), function(x, y) {
    if (is_transposed(x)) {
        return(tcrossprod(t(x)))
    }

    x0 <- get_matrix2(x)
    out <- as.matrix(crossprod(x0))

    if (use_center(x)) {
        centering <- get_center(x)
        colsums <- Matrix::colSums(x0)

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

#' @export
#' @importFrom Matrix crossprod
setMethod("crossprod", c("DeferredMatrix", "ANY"), function(x, y) {
    if (is_transposed(x)) {
        return(t(x) %*% y)
    }

    out <- as.matrix(crossprod(get_matrix2(x), y))
    if (use_center(x)) {
        if (is.null(dim(y))) {
            out <- out - get_center(x) * sum(y)
        } else {
            out <- out - outer(get_center(x), Matrix::colSums(y))
        }
    }
    if (use_scale(x)) {
        out <- out / get_scale(x)
    }
    out
})

#' @export
#' @importFrom Matrix crossprod 
setMethod("crossprod", c("ANY", "DeferredMatrix"), function(x, y) {
    if (is_transposed(y)) {
        return(t(t(y) %*% x))
    }

    out <- as.matrix(crossprod(x, get_matrix2(y)))
    if (use_center(y)) {
        if (is.null(dim(x))) {
            out <- sweep(out, 2, sum(x) * get_center(y), "-", check.margin=FALSE)
        } else {
            out <- out - outer(Matrix::colSums(x), get_center(y))
        }
    }
    if (use_scale(y)) {
        out <- sweep(out, 2, get_scale(y), "/", check.margin=FALSE)
    }
    out
})

###################################
# Transposed cross-product. 

#' @export
#' @importFrom Matrix tcrossprod
setMethod("tcrossprod", c("DeferredMatrix", "missing"), function(x, y) {
    if (is_transposed(x)) {
        return(crossprod(t(x)))
    }

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

#' @export
#' @importFrom Matrix tcrossprod
setMethod("tcrossprod", c("DeferredMatrix", "ANY"), function(x, y) {
    if (is_transposed(x)) {
        if (is.null(dim(y))) stop("non-conformable arguments")
        return(t(y %*% t(x)))
    }

    if (use_scale(x)) {
        if (!is.null(dim(y))) { # vector 'y' triggers error in 'tcrossprod' anyway.
            y <- sweep(y, 2, get_scale(x), "/", check.margin=FALSE)
        }
    }

    out <- as.matrix(tcrossprod(get_matrix2(x), y))

    if (use_center(x)) {
        out <- sweep(out, 2, as.numeric(tcrossprod(get_center(x), y)), "-", check.margin=FALSE)
    }

    out
})

#' @export
#' @importFrom Matrix tcrossprod
setMethod("tcrossprod", c("ANY", "DeferredMatrix"), function(x, y) {
    if (is_transposed(y)) {
        return(x %*% t(y))
    }

    if (use_scale(y)) {
        if (is.null(dim(x))) {
            x <- x / get_scale(y)
        } else { 
            x <- sweep(x, 2, get_scale(y), "/", check.margin=FALSE)
        }
    }

    out <- as.matrix(tcrossprod(x, get_matrix2(y)))

    if (use_center(y)) {
        out <- out - as.numeric(x %*% get_center(y))
    }

    out
})
