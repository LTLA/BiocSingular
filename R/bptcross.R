#' @importFrom BiocParallel bplapply bpnworkers
bptcross_x_by_row <- function(x, y=NULL, BPPARAM, njobs=bpnjobs_by_row(x, bpnworkers(BPPARAM))) {
    per.core <- bpsplit_by_row(x, njobs)
    if (is.null(y)) {
        y <- x 
    }
    out <- bplapply(per.core, FUN=tcrossprod, y=y, BPPARAM=BPPARAM) 
    do.call(rbind, out)
}

#' @importFrom BiocParallel bplapply bpnworkers
bptcross_y_by_row <- function(x, y, BPPARAM, njobs=bpnjobs_by_row(y, bpnworkers(BPPARAM))) {
    per.core <- bpsplit_by_row(y, njobs)
    out <- bplapply(per.core, FUN=tcrossprod, x=x, BPPARAM=BPPARAM)
    do.call(cbind, out)
}

#' @importFrom BiocParallel bpmapply bpnworkers
bptcross_by_col <- function(x, y=NULL, BPPARAM, njobs=bpnjobs_by_col(x, bpnworkers(BPPARAM))) {
    left.per.core <- bpsplit_by_col(x, njobs)
    if (is.null(y)) {
        right.per.core <- vector("list", length(left.per.core))
    } else {
        right.per.core <- bpsplit_by_col(y, njobs)
    }

    # Technically this would be slightly more memory efficient
    # if we used bpiterate to immediately Reduce jobs upon completion,
    # rather than holding all results in memory and adding them together.
    out <- bpmapply(FUN=tcrossprod, x=left.per.core, y=right.per.core, BPPARAM=BPPARAM, USE.NAMES=FALSE, SIMPLIFY=FALSE)
    Reduce("+", out)    
}

#' @importFrom BiocParallel bpnworkers
bptcross <- function(x, y=NULL, BPPARAM=SerialParam()) 
# Chooses the crossproduct function to use based on which parallelization scheme
# is the most convenient, i.e., most evenly distributes jobs among the workers.
{
    ncores <- bpnworkers(BPPARAM)
    if (ncores==1L) {
        return(tcrossprod(x, y))
    }

    if (is.null(y)) {
        x <- .matrixify_by_col(x)
        pattern <- dispatcher(x, x, ncores, x.transposed=FALSE, y.transposed=TRUE)
        if (pattern$choice=="both") {
            return(bptcross_by_col(x, BPPARAM=BPPARAM, njobs=pattern$jobs))
        } else {
            return(bptcross_x_by_row(x, BPPARAM=BPPARAM, njobs=pattern$jobs))
        }

    } else {
        x <- .matrixify_by_row(x)
        if (is.null(dim(y))) { # tcrossprod fails when 'y' is a vector.
            stop("non-conformable arguments")
        }

        pattern <- dispatcher(x, y, ncores, x.transposed=FALSE, y.transposed=TRUE)
        return(
            switch(pattern$choice,
                x=bptcross_x_by_row(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs),
                y=bptcross_y_by_row(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs),
                both=bptcross_by_col(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs)
            )
        )
    }
}
