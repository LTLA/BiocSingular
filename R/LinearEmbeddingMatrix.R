#' @title LinearEmbeddingMatrix class
#'
#' @description
#' The LinearEmbeddingMatrix class stores low-dimensional embeddings from linear dimensionality reduction methods,
#' along with information about the feature loadings and factor-level metadata.
#'
#' @param sampleFactors A matrix-like object of sample embeddings, where rows are samples and columns are factors.
#' @param featureLoadings A matrix-like object of feature loadings, where rows are features and columns are factors.
#' @param factorData A \linkS4class{DataFrame} containing factor-level information, with one row per factor.
#' @param metadata An optional list of arbitrary content describing the overall experiment.
#' 
#' @details
#' The LinearEmbeddingMatrix class is a matrix-like object where rows represent samples and columns represent factors.
#' It is designed for the storage of results from linear dimensionality reduction methods.
#' Principal components analysis (PCA) is the most obvious beneficiary,
#' but this class can also be used for other techniques like factor analysis and non-negative matrix factorization.
#' 
#' The \code{sampleFactors} slot is intended to store a low-dimensional representation of the samples, 
#' such as the principal coordinates from PCA.
#' Each row corresponds to a sample while each column corresponds to a factor.
#' 
#' The feature loadings contributing to each factor are stored in \code{featureLoadings}.
#' Each row corresponds to an input feature while each column corresponds to a factor.
#' Thus, \code{featureLoadings} should have the same number of columns as \code{sampleFactors}.
#'
#' The \code{factorData} stores additional factor-level information such as the percentage of variance explained by each factor.
#' Each row corresponds to a factor while each column represents a different metadata field.
#' Thus, it should have the same number of rows as \code{sampleFactors}.
#' 
#' This class ensures that all information related to a linear dimensionality reduction step is retained in a single object.
#' For example, feature loadings remain attached to the embedding, allowing it to be used in downstream analyses.
#'
#' @return 
#' A LinearEmbeddingMatrix object is returned from the constructor.
#' 
#' @section Getters:
#' In the following code snippets, \code{x} is a LinearEmbeddingMatrix instance.
#' \describe{
#' \item{\code{sampleFactors(x, withDimnames=TRUE)}:}{
#' Return a matrix (or matrix-like object) containing sample embeddings.
#' If \code{withDimnames=TRUE}, the dimension names of the output are the same as \code{dimnames(x)}.
#' }
#' \item{\code{factorData(x)}:}{
#' Return a \linkS4class{DataFrame} of factor-level information.
#' }
#' \item{\code{featureLoadings(x, withDimnames=TRUE)}:}{
#' Return a matrix (or matrix-like object) containing feature loadings.
#' If \code{withDimnames=TRUE}, the column names of the output are the same as \code{colnames(x)}.
#' }
#' \item{\code{dim(x)}:}{
#' Return an integer vector of length 2 containing the number of rows and columns in \code{x}.
#' }
#' \item{\code{dimnames(x)}:}{
#' Return a list of length 2 containing the row and column names.
#' Elements are either character vectors or \code{NULL}.
#' }
#' \item{\code{x$name}:}{
#' Return the field in \code{factorData(x)} with name \code{name}.
#' }
#' }
#'
#' @section Setters:
#' In the following code snippets, \code{x} is a LinearEmbeddingMatrix instance.
#' \describe{
#' \item{\code{sampleFactors(x) <- value}:}{
#' Replace the sample embeddings in \code{x} with \code{value}, 
#' a matrix (or matrix-like object) with dimensions equal to \code{dim(x)}.
#' }
#' \item{\code{factorData(x)}:}{
#' Replace the factor-level information in \code{x} with a \linkS4class{DataFrame} \code{value} 
#' that has number of rows equal to \code{nrow(x)}.
#' }
#' \item{\code{featureLoadings(x) <- value}:}{
#' Replace the feature loadings in \code{x} with \code{value},
#' a matrix (or matrix-like object) with number of columns equal to \code{ncol(x)}.
#' }
#' \item{\code{dimnames(x) <- value}:}{
#' Replace the dimension names of \code{x} with \code{value}, 
#' a list of length 2 containing the row and column names.
#' }
#' \item{\code{x$name <- value}:}{
#' Replace the field in \code{factorData(x)} with name \code{name} with \code{value}.
#' }
#' }
#'
#' @section Combining objects:
#' In the following code snippets, \code{...} contains one or more LinearEmbeddingMatrix objects.
#' \describe{
#' \item{\code{rbind(..., deparse.level=1)}:}{Returns a LinearEmbedddingMatrix 
#' where all objects in \code{...} are combined row-wise,
#' i.e., rows in successive objects are appended to the first object.
#' Refer to \code{?\link[base]{rbind}} for the interpretation of \code{deparse.level}.
#'
#' All objects in \code{...} must have the exact same values for \code{featureLoadings} and \code{factorData}.
#' }
#' \item{\code{cbind(..., deparse.level=1)}:}{Returns a LinearEmbeddingMatrix
#' where all objects in \code{...} are combined column-wise, 
#' i.e., columns in successive objects are appended to the first object.
#' Refer to \code{?\link[base]{cbind}} for the interpretation of \code{deparse.level}.
#' }
#' }
#'
#' @section Subsetting:
#' In the following code snippets, \code{x} is a LinearEmbeddingMatrix object.
#' \describe{
#' \item{\code{x[i, j, drop=TRUE]}:}{Returns a LinearEmbeddingMatrix containing the 
#' specified rows \code{i} and columns \code{j}.
#' 
#' \code{i} and \code{j} can be a logical, integer or character vector of subscripts, 
#' indicating the rows and columns respectively to retain.
#' Either can be missing, in which case subsetting is only performed in the specified dimension.
#' If both are missing, no subsetting is performed.
#' 
#' If \code{drop=TRUE} and the subsetting would produce dimensions of length 1, 
#' those dimensions are dropped and a vector is returned directly from \code{sampleFactors}.
#' }
#' \item{\code{x[i, j] <- value}:}{Replaces all data for rows \code{i} and columns {j} 
#' with the corresponding fields in a LinearEmbeddingMatrix \code{value}.
#'
#' \code{i} and \code{j} can be a logical, integer or character vector of subscripts, 
#' indicating the rows and columns respectively to replace.
#' Either can be missing, in which case replacement is only performed in the specified dimension.
#' If both are missing, \code{x} is replaced entirely with \code{value}.
#'
#' If \code{j} is specified, \code{value} is expected to have the same features in \code{featureLoadings} as \code{x}.}
#' }
#'
#' @section Other methods:
#' \code{as.matrix(x)} returns an ordinary (usually numeric) matrix of sample embeddings from a LinearEmbeddingMatrix \code{x}.
#' 
#' \code{show(object)} prints a message to screen describing the data stored in a LinearEmbeddingMatrix \code{object}.
#' 
#' @author
#' Aaron Lun, Davide Risso and Keegan Korthauer
#' 
#' @examples
#' lem <- LinearEmbeddingMatrix(matrix(rnorm(1000), ncol=5),
#'     matrix(runif(20000), ncol=5))
#' lem
#'
#' # Getting and setting:
#' sampleFactors(lem)
#' sampleFactors(lem) <- sampleFactors(lem) * -1
#' 
#' featureLoadings(lem)
#' featureLoadings(lem) <- featureLoadings(lem) * -1
#' 
#' factorData(lem)
#' factorData(lem)$whee <- 1
#' 
#' nrow(lem)
#' ncol(lem)
#' colnames(lem) <- LETTERS[seq_len(ncol(lem))]
#' as.matrix(lem)
#'
#' # Combining and subsetting:
#' rbind(lem, lem)
#' cbind(lem, lem)
#'
#' lem[1:10,]
#' lem[,1:5]
#' 
#' lem2 <- lem
#' lem2[1:10,] <- lem[11:20,]
#' 
#' @name LinearEmbeddingMatrix
#' @aliases
#' sampleFactors
#' featureLoadings
#' factorData
#' sampleFactors,LinearEmbeddingMatrix-method
#' featureLoadings,LinearEmbeddingMatrix-method
#' factorData,LinearEmbeddingMatrix-method
#' sampleFactors<-
#' featureLoadings<-
#' factorData<-
#' sampleFactors<-,LinearEmbeddingMatrix-method
#' featureLoadings<-,LinearEmbeddingMatrix-method
#' factorData<-,LinearEmbeddingMatrix-method
#' as.matrix,LinearEmbeddingMatrix-method
#' dim,LinearEmbeddingMatrix-method
#' dimnames,LinearEmbeddingMatrix-method
#' dimnames<-,LinearEmbeddingMatrix-method
#' dimnames<-,LinearEmbeddingMatrix,ANY-method
#' $,LinearEmbeddingMatrix-method
#' rbind,LinearEmbeddingMatrix-method
#' cbind,LinearEmbeddingMatrix-method
#' [,LinearEmbeddingMatrix,ANY-method
#' [,LinearEmbeddingMatrix,ANY,ANY-method
#' [,LinearEmbeddingMatrix,ANY,ANY,ANY-method
#' [<-,LinearEmbeddingMatrix,ANY,ANY,LinearEmbeddingMatrix-method
NULL

#' @export
setMethod("sampleFactors", "LinearEmbeddingMatrix", function(x, withDimnames=TRUE) {
    sf <- x@sampleFactors
    if (withDimnames) {
        colnames(sf) <- colnames(x)
    }
    sf
})

#' @export
setMethod("featureLoadings", "LinearEmbeddingMatrix", function(x, withDimnames=TRUE){
    fl <- x@featureLoadings
    if(withDimnames){
        colnames(fl) <- colnames(x)
    }
    fl
})

#' @export
setMethod("factorData", "LinearEmbeddingMatrix", function(x){
    x@factorData
})

#' @export
setReplaceMethod("sampleFactors", "LinearEmbeddingMatrix", function(x, value) {
    x@sampleFactors <- value
    validObject(x)
    x
})

#' @export
setReplaceMethod("featureLoadings", "LinearEmbeddingMatrix", function(x, value) {
    x@featureLoadings <- value
    validObject(x)
    x
})

#' @export
setReplaceMethod("factorData", "LinearEmbeddingMatrix", function(x, value) {
    x@factorData <- value
    validObject(x)
    x
})

#' @export
setMethod("$", "LinearEmbeddingMatrix", function(x, name) {
    factorData(x)[[name]]
})

#' @export
setReplaceMethod("$", "LinearEmbeddingMatrix", function(x, name, value) {
    factorData(x)[[name]] <- value
    x
})

#############################################
# Define matrix methods.

#' @export
setMethod("dim", "LinearEmbeddingMatrix", function(x) {
    dim(sampleFactors(x, withDimnames=FALSE))
})

#' @export
setMethod("dimnames", "LinearEmbeddingMatrix", function(x) {
    list(rownames(sampleFactors(x, withDimnames=FALSE)), rownames(factorData(x)))
})

#' @export
setReplaceMethod("dimnames", "LinearEmbeddingMatrix", function(x, value) {
    fd <- factorData(x)
    rownames(fd) <- value[[2]]
    factorData(x) <- fd

    sf <- sampleFactors(x, withDimnames=FALSE)
    rownames(sf) <- value[[1]]
    sampleFactors(x) <- sf
    x
})

#' @export
#' @method as.matrix LinearEmbeddingMatrix
as.matrix.LinearEmbeddingMatrix <- function(x, ...) {
    y <- sampleFactors(x)
    while (!is.matrix(y)) { # in case of sparsity or other madness.
        y <- as.matrix(y)
    }
    y
}

#' @export
setMethod("as.matrix", "LinearEmbeddingMatrix", as.matrix.LinearEmbeddingMatrix)

#############################################
# Sets the validity checker.

.le_validity <- function(object) {
    msg <- NULL

    # Check dimensions
    sf <- sampleFactors(object, withDimnames=FALSE)
    fl <- featureLoadings(object, withDimnames=FALSE)
    fd <- factorData(object, withDimnames=FALSE)

    if(length(dim(sf))!=2L || length(dim(fl))!=2L || NCOL(sf) != NCOL(fl)) {
        msg <- c(msg, "'sampleFactors' and 'featureLoadings' must have the same number of columns")
    }

    if(NROW(fd) != NCOL(sf)) {
        msg <- c(msg, "'factorData' must have one row per factor")
    }

    if (length(msg)) { return(msg) }
    return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("LinearEmbeddingMatrix", .le_validity)

#############################################
# Sets the show method.

#' @importFrom S4Vectors metadata
.le_show <- function(object) {

    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    ## metadata
    expt <- names(metadata(object))
    if (is.null(expt)) {
        expt <- character(length(metadata(object)))
    }
    scat("metadata(%d): %s\n", expt)

    ## rownames
    rownames <- rownames(object)
    if(is.null(rownames)) {
        cat("rownames: NULL\n")
    } else {
        scat("rownames(%d): %s\n", rownames(object))
    }

    ## colnames
    colnames <- colnames(object)
    if(is.null(colnames)) {
        cat("colnames: NULL\n")
    } else {
        scat("colnames(%d): %s\n", colnames(object))
    }

    ## factorData
    scat("factorData names(%d): %s\n", names(factorData(object)))
}

scat <- function(fmt, vals=character(), exdent=2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

#' @export
setMethod("show", "LinearEmbeddingMatrix", .le_show)

#############################################
# Defines a constructor

#' @export
#' @importClassesFrom S4Vectors DataFrame
#' @rdname LinearEmbeddingMatrix
LinearEmbeddingMatrix <- function(sampleFactors = matrix(nrow = 0, ncol = 0),
    featureLoadings = matrix(nrow = 0, ncol = 0),
    factorData = NULL, metadata = list()) 
{
    if (is.null(factorData)) {
        factorData <- new("DFrame", nrows = ncol(sampleFactors))
    }
    new("LinearEmbeddingMatrix",
        sampleFactors = sampleFactors,
        featureLoadings = featureLoadings,
        factorData = factorData,
        metadata = as.list(metadata))
}

#############################################
# Define subsetting methods.

.convert_subset_index <- function(subset, names) {
    if (is.character(subset)) {
        subset <- match(subset, names)
        if (any(is.na(subset))) {
            stop("subscript out of bounds")
        }
    }
    subset
}

#' @export
setMethod("[", c("LinearEmbeddingMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    temp_sf <- sampleFactors(x, withDimnames=FALSE) 
    temp_fl <- featureLoadings(x, withDimnames=FALSE)
    temp_fd <- factorData(x)

    if(!missing(i)) {
        i <- .convert_subset_index(i, rownames(x))
        temp_sf <- temp_sf[i,,drop=FALSE]
    }

    if(!missing(j)) {
        j <- .convert_subset_index(j, colnames(x))
        temp_sf <- temp_sf[,j,drop=FALSE]
        temp_fl <- temp_fl[,j,drop=FALSE]
        temp_fd <- temp_fd[j,,drop=FALSE]
    }

    # Returning a vector, a la drop=TRUE for a full matrix.
    if (any(dim(temp_sf)==1L) && drop) {
        colnames(temp_sf) <- rownames(temp_fd)
        return(drop(temp_sf))
    }

    initialize(x, sampleFactors = temp_sf, featureLoadings = temp_fl, factorData = temp_fd)
})

#' @export
setMethod("[<-", c("LinearEmbeddingMatrix", "ANY", "ANY", "LinearEmbeddingMatrix"), function(x, i, j, ..., value) {
    temp_sf <- sampleFactors(x, withDimnames=FALSE)
    temp_fl <- featureLoadings(x, withDimnames=FALSE)
    temp_fd <- factorData(x, withDimnames=FALSE)

    # i is samples
    # j is factors
    if (!missing(i)) {
        i <- .convert_subset_index(i, rownames(x))
    }
    if (!missing(j)) {
        j <- .convert_subset_index(j, colnames(x))
    }

    # Inserting sample factors.
    if (missing(i) && missing(j)) {
        temp_sf <- sampleFactors(value, withDimnames=FALSE)
    } else if (missing(i)) {
        temp_sf[,j] <- sampleFactors(value, withDimnames=FALSE)
    } else if (missing(j)) {
        temp_sf[i,] <- sampleFactors(value, withDimnames=FALSE)
    } else {
        temp_sf[i,j] <- sampleFactors(value, withDimnames=FALSE)
    }

    # Dealing with the factorData, featureLoadings.
    if (missing(i) && missing(j)) {
        temp_fl <- featureLoadings(value, withDimnames=FALSE)
        temp_fd <- factorData(value)
    } else if (!missing(j)) {
        temp_fl[,j] <- featureLoadings(value, withDimnames=FALSE)
        temp_fd[j,] <- factorData(value)
    }
   
    initialize(x, sampleFactors = temp_sf, featureLoadings = temp_fl, factorData = temp_fd)
})

#############################################
# Defining the combining methods.

#' @importFrom BiocGenerics rbind
#' @importFrom S4Vectors metadata metadata<-
setMethod("rbind", "LinearEmbeddingMatrix", function(..., deparse.level=1) {
    args <- list(...)
    x <- args[[1]]
    all_sf <- lapply(args, sampleFactors)
    all_sf <- do.call(rbind, all_sf)

    all_fd <- lapply(args, factorData)
    if (length(unique(all_fd))!=1L) {
        stop("'factorData' are not identical in 'cbind(<", class(x), ">)'")
    }
    all_fl <- lapply(args, featureLoadings)
    if (length(unique(all_fl))!=1L) {
        stop("'featureLoadings' are not identical in 'cbind(<", class(x), ">)'")
    }

    out <- initialize(x, sampleFactors=all_sf)
    metadata(out) <- unlist(lapply(args, metadata), recursive=FALSE)
    out
})

#' @importFrom BiocGenerics cbind
#' @importFrom S4Vectors metadata metadata<-
setMethod("cbind", "LinearEmbeddingMatrix", function(..., deparse.level=1) {
    args <- list(...)
    x <- args[[1]]

    all_sf <- lapply(args, sampleFactors, withDimnames=FALSE)
    all_fl <- lapply(args, featureLoadings, withDimnames=FALSE)
    all_fd <- lapply(args, factorData)

    out <- initialize(x, sampleFactors = do.call(cbind, all_sf),
        featureLoadings = do.call(cbind, all_fl),
        factorData = do.call(rbind, all_fd))
    metadata(out) <- unlist(lapply(args, metadata), recursive=FALSE)
    out
})

