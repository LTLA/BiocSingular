#' @importFrom BiocParallel bpnworkers bplapply SerialParam
bpcross_x_by_col <- function(x, y=NULL, BPPARAM=SerialParam(), njobs=NULL) {
    per.core <- bpsplit_by_col(x, bpnworkers(BPPARAM), njobs=njobs)
    if (is.null(y) && length(per.core) > 1L) { 
        y <- x 
    }
    out <- bplapply(per.core, FUN=crossprod, y=y, BPPARAM=BPPARAM) 
    do.call(rbind, out)
}

#' @importFrom BiocParallel bpnworkers bplapply SerialParam
bpcross_y_by_col <- function(x, y, BPPARAM=SerialParam(), njobs=NULL) {
    per.core <- bpsplit_by_col(y, bpnworkers(BPPARAM), njobs=njobs)
    out <- bplapply(per.core, FUN=crossprod, x=x, BPPARAM=BPPARAM)
    do.call(cbind, out)
}

#' @importFrom BiocParallel bpnworkers bpmapply SerialParam
bpcross_by_row <- function(x, y=NULL, BPPARAM=SerialParam(), njobs=NULL) {
    ncores <- bpnworkers(BPPARAM)
    if (is.null(njobs)) {
        njobs <- bpnjobs_by_row(x, ncores)
    }
    
    left.per.core <- bpsplit_by_row(x, ncores, njobs=njobs)
    if (is.null(y)) {
        right.per.core <- vector("list", ncores)
    } else {
        right.per.core <- bpsplit_by_row(y, ncores, njobs=njobs)
    }

    # Technically this would be slightly more memory efficient
    # if we used bpiterate to immediately Reduce jobs upon completion,
    # rather than holding all results in memory and adding them together.
    out <- bpmapply(FUN=crossprod, x=left.per.core, y=right.per.core, BPPARAM=BPPARAM, USE.NAMES=FALSE, SIMPLIFY=FALSE)
    Reduce("+", out)    
}

#' @importFrom BiocParallel bpnworkers
bpcross <- function(x, y=NULL, BPPARAM=SerialParam()) 
# Chooses the crossproduct function to use based on which parallelization scheme
# is the most convenient, i.e., most evenly distributes jobs among the workers.
{
    ncores <- bpnworkers(BPPARAM)
    if (ncores==1L) {
        return(crossprod(x, y))
    }

    x <- .matrixify_by_col(x)  
    njobs_by_row <- bpnjobs_by_row(x, ncores)
    njobs_by_col <- bpnjobs_by_col(x, ncores)
    row.dev <- .deviation(njobs_by_row)
    col.dev <- .deviation(njobs_by_col)

    if (is.null(y)) {
        if (row.dev < col.dev) {
            return(bpcross_by_row(x, BPPARAM=BPPARAM, njobs=njobs_by_row))
        } else {
            return(bpcross_x_by_col(x, BPPARAM=BPPARAM, njobs=njobs_by_col))
        }
    } else {
        y <- .matrixify_by_col(y)

        njobs_by_col2 <- bpnjobs_by_col(y, ncores)
        col.dev2 <- .deviation(njobs_by_col2)
        if (col.dev2 < col.dev && col.dev2 < row.dev) {
            return(bpcross_y_by_col(x, y, BPPARAM=BPPARAM, njobs=njobs_by_col2))
        } 

        if (row.dev < col.dev) {
            njobs_by_row2 <- bpnjobs_by_row(y, ncores)
            if (identical(njobs_by_row2, njobs_by_row)) { # only possible if the rows match up.
                return(bpcross_by_row(x, y, BPPARAM=BPPARAM, njobs=njobs_by_row))
            }
        }

        return(bpcross_x_by_col(x, y, BPPARAM=BPPARAM, njobs=njobs_by_col))
    }
}
