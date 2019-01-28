#' @export
bsfold <- function(object) object@fold

#' @export
bsdeferred <- function(object) object@deferred

bsargs <- function(object) object@args

#' @importFrom S4Vectors setValidity2
setValidity2("BiocSingularParam", function(object) {
    msg <- character(0)

    if (length(bsfold(object))!=1L) {
        msg <- c(msg, "'fold' should be a numeric scalar")
    }
    if (bsfold(object) < 1) {
        msg <- c(msg, "'fold' should be no less than 1")
    }

    if (length(bsdeferred(object))!=1L) {
        msg <- c(msg, "'deferred' should be a numeric scalar")
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
    cat(sprintf("cross-product fold-threshold: %.2f\n", bsfold(object)))
    cat(sprintf("deferred centering/scaling: %s\n", ifelse(bsdeferred(object), "on", "off")))
})

#' @export
#' @importFrom methods new
ExactParam <- function(deferred=FALSE, fold=5) {
    new("ExactParam", deferred=as.logical(deferred), fold=as.numeric(fold))
}

#' @export
#' @importFrom methods new
IrlbaParam <- function(deferred=FALSE, fold=5, extra.work=7, ...) {
    new("IrlbaParam", deferred=as.logical(deferred), fold=as.numeric(fold), extra.work=as.integer(extra.work), args=list(...))
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
    extra.names <- names(bsargs(object))
    if (length(extra.names) > 3) extra.names <- c(extra.names[seq_len(3)], "...")
    cat(sprintf("additional arguments(%i): %s\n", length(bsargs(object)), paste(extra.names, collapse=", ")))
})

#' @export
#' @importFrom methods new
RandomParam <- function(deferred=FALSE, fold=5, ...) {
    new("RandomParam", deferred=as.logical(deferred), fold=as.numeric(fold), args=list(...))
}

#' @export
#' @importFrom methods show 
setMethod("show", "RandomParam", function(object) {
    callNextMethod()
    extra.names <- names(bsargs(object))
    if (length(extra.names) > 3) extra.names <- c(extra.names[seq_len(3)], "...")
    cat(sprintf("additional arguments(%i): %s\n", length(bsargs(object)), paste(extra.names, collapse=", ")))
})
