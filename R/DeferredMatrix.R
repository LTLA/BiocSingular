# BiocSingular matrix class with support for deferred centering and scaling.
# Specifically, we multiply first and apply the centering/scaling afterwards.

###################################
###################################
###################################
# Constructing the seed.

#' @export
#' @import methods
setClass("DeferredMatrixSeed", slots=c(.matrix="ANY", center="numeric", scale="numeric", use_center="logical", use_scale="logical", transposed="logical"))

#' @export
#' @importFrom methods new
DeferredMatrixSeed <- function(x, center=NULL, scale=NULL) {
    use_center <- !is.null(center)
    use_scale <- !is.null(scale)
    new("DeferredMatrixSeed", .matrix=x, center=as.numeric(center), scale=as.numeric(scale), use_center=use_center, use_scale=use_scale, transposed=FALSE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("DeferredMatrixSeed", function(object) {
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
setMethod("show", "DeferredMatrixSeed", function(object) {
    cat(sprintf("%i x %i DeferredMatrixSeed object", nrow(object), ncol(object)),
        sprintf("representation: %s", class(get_matrix2(object))),
        sprintf("centering: %s", if (use_center(object)) "yes" else "no"),
        sprintf("scaling: %s", if (use_scale(object)) "yes" else "no"),
    sep="\n")
})

###################################
# Internal getters. 

get_matrix2 <- function(x) x@.matrix

get_center <- function(x) x@center

get_scale <- function(x) x@scale

use_center <- function(x) x@use_center

use_scale <- function(x) x@use_scale

is_transposed <- function(x) x@transposed

###################################
# DelayedArray support utilities. 

#' @export
setMethod("dim", "DeferredMatrixSeed", function(x) {
    d <- dim(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
setMethod("dimnames", "DeferredMatrixSeed", function(x) {
    d <- dimnames(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
#' @importFrom DelayedArray extract_array
setMethod("extract_array", "DeferredMatrixSeed", function(x, index) {
    x2 <- subset_DeferredMatrixSeed(x, index[[1]], index[[2]])        
    realize_DeferredMatrixSeed(x2)
})

###################################
# Other utilities. 

rename_DeferredMatrixSeed <- function(x, value) {
    if (is_transposed(x)) value <- rev(value)
    dimnames(x@.matrix) <- value
    x
}

#' @importFrom BiocGenerics t
transpose_DeferredMatrixSeed <- function(x) {
    x@transposed <- !is_transposed(x)
    x
}

#' @importFrom BiocGenerics t
#' @importFrom methods is
realize_DeferredMatrixSeed <- function(x, ...) {
    out <- get_matrix2(x)

    if (use_scale(x) || use_center(x)) {
        if (is(out, "DeferredMatrix")) {
            # Any '-' and '/' would collapse this to a DelayedArray, 
            # which would then call extract_array, which would then 
            # call realize_DeferredMatrixSeed, forming an infinite loop.
            # So we might as well realize it now.
            out <- realize_DeferredMatrixSeed(seed(out))
        }

        out <- t(out)
        if (nrow(out) && ncol(out)) { # avoid undefined <Matrix> - numeric(0) with Matrix classes.
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
}

#' @importFrom BiocGenerics t
subset_DeferredMatrixSeed <- function(x, i, j) {
    if (is_transposed(x)) {
        x2 <- transpose_DeferredMatrixSeed(x)
        x2 <- subset_DeferredMatrixSeed(x2, i=j, j=i)
        return(transpose_DeferredMatrixSeed(x2))
    }

    if (!is.null(i)) {
        x@.matrix <- get_matrix2(x)[i,,drop=FALSE]
    }
    
    if (!is.null(j)) {
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

    return(x)
}

###################################
# Hacks to overcome the deficiencies of "Matrix".

#' @importFrom BiocGenerics colSums
.safe_colSums <- function(x) {
    if (is(x, "Matrix")) {
        Matrix::colSums(x)
    } else {
        colSums(x)
    }    
}

#' @importFrom BiocGenerics rowSums
.safe_rowSums <- function(x) {
    if (is(x, "Matrix")) {
        Matrix::rowSums(x)
    } else {
        rowSums(x)
    }
}

###################################
###################################
###################################
# Constructing the matrix.

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
#' @import methods
setClass("DeferredMatrix",
    contains="DelayedMatrix",
    representation(seed="DeferredMatrixSeed")
)

#' @export
#' @importFrom methods new is
#' @importFrom DelayedArray DelayedArray
DeferredMatrix <- function(x, center=NULL, scale=NULL) {
    if (is(x, "DeferredMatrixSeed")) {
        seed <- x
    } else {
        seed <- DeferredMatrixSeed(x, center=center, scale=scale)
    }
    DelayedArray(seed)
}

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "DeferredMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="DeferredMatrix")
)

###################################
# Overridden utilities from DelayedArray, for efficiency.

#' @export
#' @importFrom DelayedArray DelayedArray seed
setReplaceMethod("dimnames", "DeferredMatrix", function(x, value) {
    DelayedArray(rename_DeferredMatrixSeed(seed(x), value))
})

#' @export
#' @importFrom BiocGenerics t
#' @importFrom DelayedArray DelayedArray seed
setMethod("t", "DeferredMatrix", function(x) {
    DelayedArray(transpose_DeferredMatrixSeed(seed(x)))
})

#' @export
#' @importFrom DelayedArray DelayedArray seed
setMethod("[", "DeferredMatrix", function(x, i, j, ..., drop=TRUE) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    out <- DelayedArray(subset_DeferredMatrixSeed(seed(x), i=i, j=j))

    if (drop && any(dim(out)==1L)) {
        return(drop(out))
    }
    out
})

###################################
# Basic matrix stats.

#' @export
#' @importFrom BiocGenerics colSums
#' @importFrom Matrix drop
setMethod("colSums", "DeferredMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(rowSums(t(x)))
    }

    out <- rep(1, nrow(x)) %*% x
    out <- drop(out)
    names(out) <- colnames(x)
    out
})

#' @export
#' @importFrom BiocGenerics rowSums
#' @importFrom Matrix drop
setMethod("rowSums", "DeferredMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(colSums(t(x)))
    }

    out <- x %*% rep(1, ncol(x))
    out <- drop(out)
    names(out) <- rownames(x)
    out
})

#' @export
#' @importFrom BiocGenerics colMeans
setMethod("colMeans", "DeferredMatrix", function(x, na.rm = FALSE, dims = 1L) colSums(x)/nrow(x))

#' @export
#' @importFrom BiocGenerics rowMeans
setMethod("rowMeans", "DeferredMatrix", function(x, na.rm = FALSE, dims = 1L) rowSums(x)/ncol(x))

###################################
# Matrix multiplication.

# We attempt to use operators defined for '.matrix' in the 'DeferredMatrixSeed'.
# This avoids expensive modifications such as loss of sparsity.
# Centering and scaling are factored out into separate operations.
#
# We assume that the non-'DeferredMatrix' argument is small and can be modified cheaply.
# We also assume that the matrix product is small and can be modified cheaply.
# This allows centering and scaling to be applied *after* multiplication.

#' @export
#' @importFrom BiocGenerics t
#' @importFrom DelayedArray seed DelayedArray
setMethod("%*%", c("DeferredMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        return(t(t(y) %*% t(x)))    
    }

    if (use_scale(x_seed)) {
        y <- y / get_scale(x_seed)
    }

    out <- as.matrix(get_matrix2(x_seed) %*% y)

    if (use_center(x_seed)) {
        out <- sweep(out, 2, as.numeric(get_center(x_seed) %*% y), "-", check.margin=FALSE)
    }

    DelayedArray(out)
})

#' @export
#' @importFrom BiocGenerics t
#' @importFrom DelayedArray seed DelayedArray
setMethod("%*%", c("ANY", "DeferredMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        if (!is.null(dim(x))) x <- t(x) # as vectors don't quite behave as 1-column matrices here.
        return(t(t(y) %*% x))    
    }

    out <- as.matrix(x %*% get_matrix2(y_seed))

    if (use_center(y_seed)) {
        if (is.null(dim(x))) {
            out <- out - get_center(y_seed) * sum(x)
        } else {
            out <- out - outer(.safe_rowSums(x), get_center(y_seed), "*")
        }
    }

    if (use_scale(y_seed)) {
        out <- sweep(out, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    DelayedArray(out)
})

#' @export
#' @importFrom BiocGenerics t
#' @importFrom DelayedArray seed DelayedArray
setMethod("%*%", c("DeferredMatrix", "DeferredMatrix"), function(x, y) {
    x_seed <- seed(x)
    y_seed <- seed(y)

    if (!is_transposed(x_seed)) {
        if (!is_transposed(y_seed)) {
            res <- .multiply_u2u(x_seed, y_seed)
        } else {
            res <- .multiply_u2t(x_seed, y_seed)
        }
    } else {
        if (!is_transposed(y_seed)) {
            res <- .multiply_t2u(x_seed, y_seed)
        } else {
            res <- .multiply_u2u(y_seed, x_seed)
            res <- t(res)
        }
    }

    DelayedArray(res)
})

###################################
# DefMat %*% DefMat utilities.

#' @importFrom Matrix drop
.multiply_u2u <- function(x_seed, y_seed) 
# Considering the problem of (X - C_x)S_x (Y - C_y)S_y.
{
    # Computing X S_x Y S_y
    x0 <- get_matrix2(x_seed)
    if (use_scale(x_seed)) {
        x0 <- DeferredMatrix(x0, scale=get_scale(x_seed))
    } 

    result <- as.matrix(x0 %*% get_matrix2(y_seed))
    if (use_scale(y_seed)) {
        result <- sweep(result, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    # Computing C_x S_x Y S_y, and subtracting it from 'result'.
    if (use_center(x_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }

        component2 <- drop(x.center %*% get_matrix2(y_seed))
        if (use_scale(y_seed)) {
            component2 <- component2 / get_scale(y_seed)
        }

        result <- sweep(result, 2, component2, "-", check.margin=FALSE)
    }

    # Computing C_x S_x C_y S_y, and adding it to 'result'.
    if (use_center(x_seed) && use_center(y_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }

        y.center <- get_center(y_seed)
        if (use_scale(y_seed)) {
            y.center <- y.center / get_scale(y_seed)
        }

        component4 <- sum(x.center) * y.center
        result <- sweep(result, 2, component4, "+", check.margin=FALSE)
    }

    # Computing X S_x C_y S_y, and subtracting it from 'result'.
    # This is done last to avoid subtracting large values.
    if (use_center(y_seed)) {
        y.center <- get_center(y_seed)
        if (use_scale(y_seed)) {
            y.center <- y.center / get_scale(y_seed)
        }

        component3 <- outer(.safe_rowSums(x0), y.center)
        result <- result - component3
    }

    result
}

#' @importFrom Matrix tcrossprod drop
.multiply_u2t <- function(x_seed, y_seed) 
# Considering the problem of (X - C_x)S_x S_y(Y' - C_y')
{
    # Computing X S_x S_y Y'
    x0 <- get_matrix2(x_seed)
    if (use_scale(x_seed) || use_scale(y_seed)) {
        scaling <- 1
        if (use_scale(x_seed)) {
            scaling <- scaling * get_scale(x_seed)
        }
        if (use_scale(y_seed)) {
            scaling <- scaling * get_scale(y_seed)
        }
        x0 <- DeferredMatrix(x0, scale=scaling)
    }
    result <- as.matrix(tcrossprod(x0, get_matrix2(y_seed)))

    # Computing C_x S_x S_y Y', and subtracting it from 'result'.
    if (use_center(x_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }
        if (use_scale(y_seed)) {
            x.center <- x.center / get_scale(y_seed)
        }

        component2 <- drop(tcrossprod(x.center, get_matrix2(y_seed)))
        result <- sweep(result, 2, component2, "-", check.margin=FALSE)
    }

    # Computing C_x S_x S_y C_y', and adding it to 'result'.
    if (use_center(x_seed) && use_center(y_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }

        y.center <- get_center(y_seed)
        if (use_scale(y_seed)) {
            y.center <- y.center / get_scale(y_seed)
        }

        component4 <- sum(x.center*y.center)
        result <- result + component4
    }

    # Computing X S_x S_y C_y', and subtracting it from 'result'.
    # This is done last to avoid subtracting large values.
    if (use_center(y_seed)) {
        component3 <- drop(x0 %*% get_center(y_seed))
        result <- result - component3
    }

    result
}

#' @importFrom Matrix crossprod
.multiply_t2u <- function(x_seed, y_seed) 
# Considering the problem of S_x(X' - C_x') (Y - C_y)S_y
{
    # C mputing X' Y 
    result <- as.matrix(crossprod(get_matrix2(x_seed), get_matrix2(y_seed)))

    # Computing C_x' Y, and subtracting it from 'result'.
    if (use_center(x_seed)) {
        x.center <- get_center(x_seed)
        component2 <- outer(x.center, .safe_colSums(get_matrix2(y_seed)))
        result <- result - component2
    }

    # Computing C_x' C_y, and adding it to 'result'.
    if (use_center(x_seed) && use_center(y_seed)) {
        x.center <- get_center(x_seed)
        y.center <- get_center(y_seed)
        component4 <- outer(x.center, y.center) * nrow(y_seed)
        result <- result + component4
    }

    # Computing X' C_y, and subtracting it from 'result'.
    # This is done last to avoid subtracting large values.
    if (use_center(y_seed)) {
        component3 <- outer(.safe_colSums(get_matrix2(x_seed)), get_center(y_seed))
        result <- result - component3
    }

    if (use_scale(x_seed)) {
        result <- result / get_scale(x_seed)
    } 
    if (use_scale(y_seed)) {
        result <- sweep(result, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    result
}

###################################
# Cross-product. 

#' @export
#' @importFrom Matrix crossprod tcrossprod 
#' @importFrom DelayedArray seed DelayedArray
setMethod("crossprod", c("DeferredMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        return(tcrossprod(t(x)))
    }

    x0 <- get_matrix2(x_seed)
    out <- as.matrix(crossprod(x0))

    if (use_center(x_seed)) {
        centering <- get_center(x_seed)
        colsums <- .safe_colSums(x0)

        # Minus, then add, then minus, to mitigate cancellation.
        out <- out - outer(centering, colsums)
        out <- out + outer(centering, centering) * nrow(x0)
        out <- out - outer(colsums, centering)
    }

    if (use_scale(x_seed)) {
        scaling <- get_scale(x_seed)
        out <- sweep(out / scaling, 2, scaling, "/", check.margin=FALSE)
    }

    DelayedArray(out)
})

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray seed DelayedArray
#' @importFrom BiocGenerics t
setMethod("crossprod", c("DeferredMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        return(t(x) %*% y)
    }

    out <- as.matrix(crossprod(get_matrix2(x_seed), y))

    if (use_center(x_seed)) {
        if (is.null(dim(y))) {
            out <- out - get_center(x_seed) * sum(y)
        } else {
            out <- out - outer(get_center(x_seed), .safe_colSums(y))
        }
    }
    
    if (use_scale(x_seed)) {
        out <- out / get_scale(x_seed)
    }

    DelayedArray(out)
})

#' @export
#' @importFrom Matrix crossprod 
#' @importFrom DelayedArray seed DelayedArray
#' @importFrom BiocGenerics t
setMethod("crossprod", c("ANY", "DeferredMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        return(t(t(y) %*% x))
    }

    out <- as.matrix(crossprod(x, get_matrix2(y_seed)))

    if (use_center(y_seed)) {
        if (is.null(dim(x))) {
            out <- sweep(out, 2, sum(x) * get_center(y_seed), "-", check.margin=FALSE)
        } else {
            out <- out - outer(.safe_colSums(x), get_center(y_seed))
        }
    }

    if (use_scale(y_seed)) {
        out <- sweep(out, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    DelayedArray(out)
})

#' @export
setMethod("crossprod", c("DeferredMatrix", "DeferredMatrix"), function(x, y) {
    t(x) %*% y
})

###################################
# Transposed cross-product. 

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray seed DelayedArray
#' @importFrom BiocGenerics t
setMethod("tcrossprod", c("DeferredMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        return(crossprod(t(x)))
    }

    x0 <- get_matrix2(x_seed)

    if (use_scale(x_seed)) {
        out <- as.matrix(.internal_tcrossprod(x0, get_scale(x_seed)))
    } else {
        out <- as.matrix(tcrossprod(x0))
    }
    
    if (use_center(x_seed)) {
        centering <- get_center(x_seed)

        if (use_scale(x_seed)) {
            centering <- centering / get_scale(x_seed)
            extra <- centering / get_scale(x_seed)
        } else {
            extra <- centering
        }
            
        # With scaling, the use of 'extra' mimics sweep(x0, 2, get_scale(x), "/"),
        # except that the scaling is applied to 'centering' rather than directly to 'x0'.
        # Without scaling, 'extra' and 'centering' are interchangeable.
        component <- tcrossprod(extra, x0)

        # Minus, then add, then minus, to mitigate cancellation.
        out <- sweep(out, 2, as.numeric(component), "-", check.margin=FALSE)
        out <- out + sum(centering^2)
        out <- out - as.numeric(x0 %*% extra)
    }

    DelayedArray(out)
})

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray seed DelayedArray
#' @importFrom BiocGenerics t
setMethod("tcrossprod", c("DeferredMatrix", "ANY"), function(x, y) {
    if (is.null(dim(y))) { # for consistency with base::tcrossprod.
        stop("non-conformable arguments")
    }

    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        return(t(y %*% t(x)))
    }

    if (use_scale(x_seed)) {
        # 'y' cannot be a vector anymore, due to the check above.
        y <- sweep(y, 2, get_scale(x_seed), "/", check.margin=FALSE)
    }

    out <- as.matrix(tcrossprod(get_matrix2(x_seed), y))

    if (use_center(x_seed)) {
        out <- sweep(out, 2, as.numeric(tcrossprod(get_center(x_seed), y)), "-", check.margin=FALSE)
    }

    DelayedArray(out)
})

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray seed DelayedArray
#' @importFrom BiocGenerics t
setMethod("tcrossprod", c("ANY", "DeferredMatrix"), function(x, y) {
    y_seed <- seed(y) 
    if (is_transposed(y_seed)) {
        return(x %*% t(y))
    }

    if (use_scale(y_seed)) {
        if (is.null(dim(x))) {
            x <- x / get_scale(y_seed)
        } else { 
            x <- sweep(x, 2, get_scale(y_seed), "/", check.margin=FALSE)
        }
    }

    out <- as.matrix(tcrossprod(x, get_matrix2(y_seed)))

    if (use_center(y_seed)) {
        out <- out - as.numeric(x %*% get_center(y_seed))
    }

    DelayedArray(out)
})

#' @export
setMethod("tcrossprod", c("DeferredMatrix", "DeferredMatrix"), function(x, y) {
    x %*% t(y)
})

###################################
# Extra code for corner-case calculations of the transposed cross-product.

#' @importFrom DelayedArray seed DelayedArray
.update_scale <- function(x, s) {
    x_seed <- seed(x)
    if (use_scale(x_seed)) {
        s <- s * get_scale(x_seed)
    }
    x_seed@scale <- s
    x_seed@use_scale <- TRUE
    DelayedArray(x_seed)
}

#' @importFrom Matrix tcrossprod 
#' @importFrom methods is
#' @importFrom DelayedArray seed
.internal_tcrossprod <- function(x, scale.) 
# Computes tcrossprod(sweep(x, 2, scale, "/")) when 'x' is a DeferredMatrix.
# 'scale' can be assumed to be non-NULL here.
# This will always return a dense ordinary matrix.
{
    if (!is(x, "DeferredMatrix")) {
        x <- sweep(x, 2, scale., "/", check.margin=FALSE) 
        return(as.matrix(tcrossprod(x)))
    }

    x_seed <- seed(x)
    if (!is_transposed(x_seed)) {
        x <- .update_scale(x, scale.)
        return(as.matrix(tcrossprod(x)))
    }

    inner <- get_matrix2(x_seed)
    if (is(inner, "DeferredMatrix")) {
        if (is_transposed(seed(inner))) {
            component1 <- as.matrix(crossprod(.update_scale(inner, scale.)))
        } else {
            component1 <- .internal_tcrossprod(t(inner), scale.) # recurses. 
        }
    } else {
        component1 <- as.matrix(crossprod(inner/scale.))
    }
           
    if (use_center(x_seed)) {
        centering <- get_center(x_seed)
        component2 <- .internal_mult_special(centering, scale., inner)
        component3 <- t(component2)
        component4 <- outer(centering, centering) * sum(1/scale.^2)
        final <- (component1 - component2) + (component4 - component3)
    } else {
        final <- component1
    }

    if (use_scale(x_seed)) {
        x.scale <- get_scale(x_seed)
        final <- final / x.scale
        final <- sweep(final, 2, x.scale, "/", check.margin=FALSE) 
    }

    final 
}

#' @importFrom methods is
#' @importFrom DelayedArray seed
.internal_mult_special <- function(center, scale., Z)
# Computes C^T * S^2 * Z where C is a matrix of 'centers' copied byrow=TRUE;
# S is a diagonal matrix filled with '1/scale'; and 'Z' is a DeferredMatrix.
# This will always return a dense ordinary matrix.
{
    if (!is(Z, "DeferredMatrix")) {
        return(outer(center, .safe_colSums(Z / scale.^2)))
    }
    
    Z_seed <- seed(Z)
    if (is_transposed(Z_seed)) {
        Z <- .update_scale(Z, scale.^2)
        return(outer(center, .safe_colSums(Z)))
    }

    output <- .internal_mult_special(center, scale., get_matrix2(Z_seed)) # recurses.

    if (use_center(Z_seed)) {
        output <- output - outer(center, get_center(Z_seed)) * sum(1/scale.^2)
    }

    if (use_scale(Z_seed)) {
        output <- sweep(output, 2, get_scale(Z_seed), "/")
    }

    output
}
