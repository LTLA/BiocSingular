# BiocSingular matrix class for a reconstructed low-rank matrix.
# Useful for representing the low-rank matrix without actually making it.

###################################
###################################
###################################
# Constructing the seed.

#' @export
#' @importFrom methods new is
LowRankMatrixSeed <- function(rotation, components) {
    if (missing(rotation)) {
        if (missing(components)) {
            rotation <- components <- matrix(0, 0, 0)
        } else {
            rotation <- matrix(0, 0, ncol(components))
        }
    } else if (is(rotation, "LowRankMatrixSeed")) {
        return(rotation) # directly returning the seed already.
    } else if (missing(components)) {
        components <- matrix(0, 0, ncol(rotation))
    }

    new("LowRankMatrixSeed", rotation=rotation, components=components)
}

#' @importFrom S4Vectors setValidity2
setValidity2("LowRankMatrixSeed", function(object) {
    msg <- character(0)

    R <- get_rotation(object)
    C <- get_components(object)

    if (length(dim(R))!=2L || length(dim(C))!=2L) {
        msg <- c(msg, "'components' and 'rotation' must be matrix-like objects")
    } else if (ncol(R) != ncol(C)) {
        msg <- c(msg, "number of columns in 'components' and 'rotation' must be the same");
    }

    if (length(msg)) {
        return(msg)
    } 
    return(TRUE)
})

#' @export
#' @importFrom methods show
setMethod("show", "LowRankMatrixSeed", function(object) {
    cat(sprintf("%i x %i LowRankMatrixSeed object", nrow(object), ncol(object)),
    sep="\n")
})

###################################
# Internal getters.

get_rotation <- function(x) x@rotation

get_components <- function(x) x@components

###################################
# DelayedArray support utilities. 

#' @export
setMethod("dim", "LowRankMatrixSeed", function(x) {
    c(nrow(get_rotation(x)), nrow(get_components(x)))
})

#' @export
setMethod("dimnames", "LowRankMatrixSeed", function(x) {
    list(rownames(get_rotation(x)), rownames(get_components(x)))
})

#' @export
#' @importFrom DelayedArray extract_array
#' @importFrom Matrix tcrossprod
setMethod("extract_array", "LowRankMatrixSeed", function(x, index) {
	x2 <- subset_LowRankMatrixSeed(x, index[[1]], index[[2]])
    as.matrix(tcrossprod(get_rotation(x2), get_components(x2)))
})

###################################
# Additional utilities for efficiency.

subset_LowRankMatrixSeed <- function(x, i, j) {
    C <- get_components(x)
    R <- get_rotation(x)

    if (!is.null(i)) {
        R <- R[i,,drop=FALSE]
    }
    
    if (!is.null(j)) {
        C <- C[j,,drop=FALSE]
    }

    LowRankMatrixSeed(R, C)
}

transpose_LowRankMatrixSeed <- function(x) {
    C <- get_components(x)
    R <- get_rotation(x)
    LowRankMatrixSeed(C, R)
}

rename_LowRankMatrixSeed <- function(x, value) {
    R <- get_rotation(x)
    rownames(R) <- value[[1]]
    C <- get_components(x)
    rownames(C) <- value[[2]]
    LowRankMatrixSeed(R, C)
}

###################################
###################################
###################################
# Constructing the matrix.

#' @export
#' @importFrom DelayedArray DelayedArray
LowRankMatrix <- function(rotation, components) {
    DelayedArray(LowRankMatrixSeed(rotation, components))
}

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "LowRankMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="LowRankMatrix")
)

###################################
# Overridden utilities from DelayedArray, for efficiency.

#' @export
#' @importFrom DelayedArray DelayedArray seed
setReplaceMethod("dimnames", "LowRankMatrix", function(x, value) {
    DelayedArray(rename_LowRankMatrixSeed(seed(x), value))
})

#' @export
#' @importFrom BiocGenerics t
#' @importFrom DelayedArray DelayedArray seed
setMethod("t", "LowRankMatrix", function(x) {
    DelayedArray(transpose_LowRankMatrixSeed(seed(x)))
})

#' @export
#' @importFrom DelayedArray DelayedArray seed
setMethod("[", "LowRankMatrix", function(x, i, j, ..., drop=TRUE) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    out <- DelayedArray(subset_LowRankMatrixSeed(seed(x), i=i, j=j))

    if (drop && any(dim(out)==1L)) {
        return(drop(out))
    }
    out
})
