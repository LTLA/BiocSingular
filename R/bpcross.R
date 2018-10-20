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
