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
setReplaceMethod("dimnames", "DeferredMatrix", function(x, value) {
    if (is_transposed(x)) value <- rev(value)
    dimnames(x@.matrix) <- value
    x
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
#' @importFrom methods is
setMethod("as.matrix", "DeferredMatrix", function(x) {
    out <- get_matrix2(x)

    if (use_scale(x) || use_center(x)) {
        if (is(out, "DeferredMatrix")) {
            # can't define '-' and '/' in general, so we might as well realize it now.
            out <- as.matrix(out)
        }

        out <- t(out)
        if (nrow(out)) { # avoid undefined <Matrix> - numeric(0) with Matrix classes.
            if (use_center(x)) {
                out <- out - get_center(x)
            }
            if (use_scale(x)) {
                out <- out / get_scale(x)
            }
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
            if (is.character(j)) {
                j <- match(j, colnames(x))
            }

            x@.matrix <- get_matrix2(x)[,j,drop=FALSE]
            
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
# Matrix stats.

#' @export
#' @importFrom BiocGenerics colSums
setMethod("colSums", "DeferredMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(x)) {
        return(rowSums(t(x)))
    }

    out <- rep(1, nrow(x)) %*% x
    out <- drop(out)
    names(out) <- colnames(x)
    out
})

#' @export
#' @importFrom BiocGenerics rowSums
setMethod("rowSums", "DeferredMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(x)) {
        return(colSums(t(x)))
    }

    out <- x %*% rep(1, ncol(x))
    out <- drop(out)
    names(out) <- rownames(x)
    out
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
        if (is(x0, "DeferredMatrix")) { # TODO: fix when BiocGenerics and Matrix stop fighting each other.
            colsums <- colSums(x0)
        } else {
            colsums <- Matrix::colSums(x0)
        }

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

#' @export
setMethod("crossprod", c("DeferredMatrix", "DeferredMatrix"), function(x, y) {
    x %*% y # leads to error above.
})

###################################
# Transposed cross-product. 

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom methods is
setMethod("tcrossprod", c("DeferredMatrix", "missing"), function(x, y) {
    if (is_transposed(x)) {
        return(crossprod(t(x)))
    }

    new.x <- get_matrix2(x)
    if (use_scale(x)) {
        out <- as.matrix(.internal_tcrossprod(new.x, get_scale(x)))
    } else {
        out <- as.matrix(tcrossprod(new.x))
    }
    
    if (use_center(x)) {
        centering <- get_center(x)

        if (use_scale(x)) {
            centering <- centering / get_scale(x)
            extra <- centering / get_scale(x)
        } else {
            extra <- centering
        }
            
        # With scaling, the use of 'extra' mimics sweep(new.x, 2, get_scale(x), "/"),
        # except that the scaling is applied to 'centering' rather than directly to 'new.x'.
        # Without scaling, 'extra' and 'centering' are interchangeable.
        component <- tcrossprod(extra, new.x) 

        # Minus, then add, then minus, to mitigate cancellation.
        out <- sweep(out, 2, as.numeric(component), "-", check.margin=FALSE)
        out <- out + sum(centering^2)
        out <- out - as.numeric(new.x %*% extra)
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

#' @export
setMethod("tcrossprod", c("DeferredMatrix", "DeferredMatrix"), function(x, y) {
    x %*% y # leads to error above.
})

###################################
# Extra code for corner-case calculations of the transposed cross-product.

.update_scale <- function(x, s) {
    if (use_scale(x)) {
        s <- s * get_scale(x)
    }
    x@scale <- s
    x@use_scale <- TRUE
    x
}

#' @importFrom Matrix tcrossprod 
#' @importFrom methods is
.internal_tcrossprod <- function(x, scale.) 
# Tries to compute tcrossprod(sweep(x, 2, scale, "/")) when 'x' is a DeferredMatrix.
# 'scale' can be assumed to be non-NULL here.
{
    if (!is(x, "DeferredMatrix")) {
        x <- sweep(x, 2, scale., "/", check.margin=FALSE) 
        return(tcrossprod(x))
    }

    if (!is_transposed(x)) {
        x <- .update_scale(x, scale.)
        return(tcrossprod(x))
    }

    inner <- get_matrix2(x)
    if (is(inner, "DeferredMatrix")) {
        if (is_transposed(inner)) {
            component1 <- crossprod(.update_scale(inner, scale.))
        } else {
            component1 <- .internal_tcrossprod(t(inner), scale.) # recurses. 
        }
    } else {
        component1 <- crossprod(inner/scale.)
    }
           
    if (use_center(x)) {
        centering <- get_center(x)
        component2 <- .internal_mult_special(centering, scale., inner)
        component3 <- t(component2)
        component4 <- outer(centering, centering) * sum(1/scale.^2)
        final <- (component1 - component2) + (component4 - component3)
    } else {
        final <- component1
    }

    if (use_scale(x)) {
        x.scale <- get_scale(x)
        final <- final / x.scale
        final <- sweep(final, 2, x.scale, "/", check.margin=FALSE) 
    }

    final 
}

.internal_mult_special <- function(center, scale., Z)
# Computes C^T * S^2 * Z where C is a matrix of 'centers' copied byrow=TRUE;
# S is a diagonal matrix filled with 'scale'; and 'Z' is a DeferredMatrix.
{
    if (!is(Z, "DeferredMatrix")) {
        return(outer(center, colSums(Z/scale.^2)))
    }
        
    if (is_transposed(Z)) {
        Z <- .update_scale(Z, scale.^2)
        return(outer(center, colSums(Z)))
    }

    output <- .internal_mult_special(center, scale., get_matrix2(Z)) # recurses.

    if (use_center(Z)) {
        output <- output - outer(center, get_center(Z)) * sum(1/scale.^2)
    }

    if (use_scale(Z)) {
        output <- sweep(output, 2, get_scale(Z), "/")
    }

    output
}
