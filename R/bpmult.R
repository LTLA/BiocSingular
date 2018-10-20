#' @importFrom BiocParallel bpnworkers bplapply SerialParam
bpmult_x_by_row <- function(x, y, BPPARAM, njobs=NULL) {
    per.core <- bpsplit_by_row(x, bpnworkers(BPPARAM), njobs=njobs)
    out <- bplapply(per.core, FUN=.internal_mult, right=y, BPPARAM=BPPARAM) 
    do.call(rbind, out)
}

#' @importFrom BiocParallel bpnworkers bplapply SerialParam
bpmult_y_by_col <- function(x, y, BPPARAM, njobs=NULL) {
    per.core <- bpsplit_by_col(y, bpnworkers(BPPARAM), njobs=njobs)
    out <- bplapply(per.core, FUN=.internal_mult, left=x, BPPARAM=BPPARAM) 
    do.call(cbind, out)
}

#' @importFrom BiocParallel bpnworkers bpmapply SerialParam
bpmult_x_by_col <- function(x, y, BPPARAM, njobs=NULL) {
    ncores <- bpnworkers(BPPARAM)
    if (is.null(njobs)) {
        njobs <- bpnjobs_by_col(x, ncores)
    }

    left.per.core <- bpsplit_by_col(x, ncores, njobs=njobs)
    right.per.core <- bpsplit_by_row(y, ncores, njobs=njobs)
 
    # Technically this would be slightly more memory efficient
    # if we used bpiterate to immediately Reduce jobs upon completion,
    # rather than holding all results in memory and adding them together.
    out <- bpmapply(FUN=.internal_mult, left=left.per.core, right=right.per.core, BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    Reduce("+", out)    
}

#' @importFrom BiocParallel bpnworkers SerialParam
bpmult <- function(x, y, BPPARAM=SerialParam()) 
# Parallelizing the multiplication by splitting the work by rows of 'x' or columns of 'y',
# depending on which one yields a smaller per-worker serialization.
{ 
    ncores <- bpnworkers(BPPARAM)
    if (ncores==1L) {
        return(x %*% y)
    }

    x <- .matrixify_by_row(x)
    y <- .matrixify_by_col(y)

    njobs_by_row <- bpnjobs_by_row(x, ncores)
    njobs_by_col <- bpnjobs_by_col(x, ncores)
    njobs_by_col2 <- bpnjobs_by_col(y, ncores)
    
    row.dev <- .deviation(njobs_by_row)
    col.dev <- .deviation(njobs_by_col)
    col.dev2 <- .deviation(njobs_by_col2)

    if (col.dev2 < col.dev && col.dev2 < row.dev) {
        return(bpmult_y_by_col(x, y, BPPARAM=BPPARAM))
    }

    if (col.dev < row.dev) {
        njobs_by_row2 <- bpnjobs_by_row(y, ncores)
        if (identical(njobs_by_col, njobs_by_row2)) { # only possible if the columns of 'x' match up to rows of 'y'.
            return(bpmult_x_by_col(x, y, BPPARAM=BPPARAM))
        }
    }

    return(bpmult_x_by_row(x, y, BPPARAM=BPPARAM))
}

.internal_mult <- function(left, right) 
# Avoid problems with using 'y' in lapply on '%*%' directly.
{
    left %*% right
}

.matrixify_by_row <- function(x) {
    if (is.null(dim(x))) { 
        dim(x) <- c(1L, length(x))
    }
    x
}

.matrixify_by_col <- function(x) {
    if (is.null(dim(x))) { 
        dim(x) <- c(length(x), 1L)
    }
    x
}
