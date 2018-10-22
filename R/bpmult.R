#' @importFrom BiocParallel bplapply bpnworkers
bpmult_x_by_row <- function(x, y, BPPARAM, njobs=bpnjobs_by_row(x, bpnworkers(BPPARAM))) {
    per.core <- bpsplit_by_row(x, njobs)
    out <- bplapply(per.core, FUN=.internal_mult, right=y, BPPARAM=BPPARAM) 
    do.call(rbind, out)
}

#' @importFrom BiocParallel bplapply bpnworkers
bpmult_y_by_col <- function(x, y, BPPARAM, njobs=bpnjobs_by_col(y, bpnworkers(BPPARAM))) {
    per.core <- bpsplit_by_col(y, njobs)
    out <- bplapply(per.core, FUN=.internal_mult, left=x, BPPARAM=BPPARAM) 
    do.call(cbind, out)
}

#' @importFrom BiocParallel bpmapply bpnworkers
bpmult_x_by_col <- function(x, y, BPPARAM, njobs=bpnjobs_by_col(x, bpnworkers(BPPARAM))) {
    left.per.core <- bpsplit_by_col(x, njobs)
    right.per.core <- bpsplit_by_row(y, njobs)
 
    # Technically this would be slightly more memory efficient
    # if we used bpiterate to immediately Reduce jobs upon completion,
    # rather than holding all results in memory and adding them together.
    out <- bpmapply(FUN=.internal_mult, left=left.per.core, right=right.per.core, BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    Reduce("+", out)    
}

#' @importFrom BiocParallel bpnworkers SerialParam
bpmult <- function(x, y, BPPARAM=SerialParam()) 
# Choosing how to parallelize the matrix multiplication.
{ 
    ncores <- bpnworkers(BPPARAM)
    if (ncores==1L) {
        return(.internal_mult(x, y))
    }

    x <- .matrixify_by_row(x)
    y <- .matrixify_by_col(y)
    pattern <- dispatcher(x, y, ncores, x.transposed=FALSE, y.transposed=FALSE)

    switch(pattern$choice,
        x=bpmult_x_by_row(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs),
        y=bpmult_y_by_col(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs),
        both=bpmult_x_by_col(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs)
    )
}

.internal_mult <- function(left, right) 
# Avoid problems with using 'y' in lapply on '%*%' directly.
# Also ensure that the output is a dense ordinary matrix.
{
    as.matrix(left %*% right)
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
