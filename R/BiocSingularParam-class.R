bsp_fold <- function(object) object@fold
bsp_args <- function(object) object@args

setValidity("BiocSingularParam", function(object) {
    msg <- character(0)
    if (bsp_fold(object) < 1) {
        msg <- c(msg, "'fold' should be no less than 1")
    }
    if (length(msg)) {
        return(msg)
    }
    return(TRUE)
})

#' @export
#' @importFrom methods show 
setMethod("show", "BiocSingularParam", function(object) {
    cat(sprintf("class: %s\n", class(object)))
    cat(sprintf("cross-product fold-threshold: %.2f\n", bsp_fold(object)))
})

#' @export
#' @importFrom methods new
ExactParam <- function(fold=5) {
    new("ExactParam", fold=as.numeric(fold))
}

#' @export
#' @importFrom methods new
IrlbaParam <- function(fold=5, extra.work=7, ...) {
    new("IrlbaParam", fold=as.numeric(fold), extra.work=as.integer(extra.work), args=list(...))
}

ip_extra <- function(object) object@extra.work

setValidity("IrlbaParam", function(object) {
    msg <- character(0)
    if (ip_extra(object) < 0L) {
        msg <- c(msg, "'extra.work' should be non-negative")
    }
    if (length(msg)) {
        return(msg)
    }
    return(TRUE)
})

#' @export
#' @importFrom methods show 
setMethod("show", "IrlbaParam", function(object) {
    callNextMethod()
    cat(sprintf("extra workspace: %i\n", ip_extra(object)))
    extra.names <- names(bsp_args(object))
    if (length(extra.names) > 3) extra.names <- c(extra.names[seq_len(3)], "...")
    cat(sprintf("additional arguments(%i): %s\n", length(bsp_args(object)), paste(extra.names, collapse=", ")))
})

#' @export
#' @importFrom methods new
RandomParam <- function(fold=5, ...) {
    new("RandomParam", fold=as.numeric(fold), args=list(...))
}

#' @export
#' @importFrom methods show 
setMethod("show", "RandomParam", function(object) {
    callNextMethod()
    extra.names <- names(bsp_args(object))
    if (length(extra.names) > 3) extra.names <- c(extra.names[seq_len(3)], "...")
    cat(sprintf("additional arguments(%i): %s\n", length(bsp_args(object)), paste(extra.names, collapse=", ")))
})
