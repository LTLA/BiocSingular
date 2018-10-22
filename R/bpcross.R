#' @importFrom BiocParallel bplapply bpnworkers
bpcross_x_by_col <- function(x, y=NULL, BPPARAM, njobs=bpnjobs_by_col(x, bpnworkers(BPPARAM))) {
    per.core <- bpsplit_by_col(x, njobs)
    if (is.null(y)) { 
        y <- x 
    }
    out <- bplapply(per.core, FUN=.internal_crossprod, y=y, BPPARAM=BPPARAM) 
    do.call(rbind, out)
}

#' @importFrom BiocParallel bplapply bpnworkers
bpcross_y_by_col <- function(x, y, BPPARAM, njobs=bpnjobs_by_col(y, bpnworkers(BPPARAM))) {
    per.core <- bpsplit_by_col(y, njobs)
    out <- bplapply(per.core, FUN=.internal_crossprod, x=x, BPPARAM=BPPARAM)
    do.call(cbind, out)
}

#' @importFrom BiocParallel bpmapply 
bpcross_by_row <- function(x, y=NULL, BPPARAM, njobs=bpnjobs_by_row(x, bpnworkers(BPPARAM))) {
    left.per.core <- bpsplit_by_row(x, njobs)
    if (is.null(y)) {
        right.per.core <- vector("list", length(left.per.core))
    } else {
        right.per.core <- bpsplit_by_row(y, njobs)
    }

    # Technically this would be slightly more memory efficient
    # if we used bpiterate to immediately Reduce jobs upon completion,
    # rather than holding all results in memory and adding them together.
    out <- bpmapply(FUN=.internal_crossprod, x=left.per.core, y=right.per.core, BPPARAM=BPPARAM, USE.NAMES=FALSE, SIMPLIFY=FALSE)
    Reduce("+", out)    
}

#' @importFrom BiocParallel bpnworkers
bpcross <- function(x, y=NULL, BPPARAM=SerialParam()) 
# Chooses the crossproduct function to use based on which parallelization scheme
# is the most convenient, i.e., most evenly distributes jobs among the workers.
{
    ncores <- bpnworkers(BPPARAM)
    if (ncores==1L) {
        return(.internal_crossprod(x, y))
    }

    x <- .matrixify_by_col(x)  

    if (is.null(y)) {
        pattern <- dispatcher(x, x, ncores, x.transposed=TRUE, y.transposed=FALSE)
        if (pattern$choice == "both") {
            return(bpcross_by_row(x, BPPARAM=BPPARAM, njobs=pattern$jobs))
        } else {
            return(bpcross_x_by_col(x, BPPARAM=BPPARAM, njobs=pattern$jobs))
        }

    } else {
        y <- .matrixify_by_col(y)
        pattern <- dispatcher(x, y, ncores, x.transposed=TRUE, y.transposed=FALSE)
        return(
            switch(pattern$choice,
                x=bpcross_x_by_col(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs),
                y=bpcross_y_by_col(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs),
                both=bpcross_by_row(x, y, BPPARAM=BPPARAM, njobs=pattern$jobs)
            )
        )
    }
}

.internal_crossprod <- function(x, y) 
# Ensure that the output is a dense ordinary matrix.
{
    as.matrix(crossprod(x, y))
}
