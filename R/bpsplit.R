bpnjobs_by_row <- function(x, ncores) {
    inc <- ceiling(nrow(x)/ncores)
    out <- rep(inc, ncores)
    out[ncores] <- nrow(x) - (ncores - 1L) * inc
    out
}

bpnjobs_by_col <- function(x, ncores) {
    inc <- ceiling(ncol(x)/ncores)
    out <- rep(inc, ncores)
    out[ncores] <- ncol(x) - (ncores - 1L) * inc
    out
}

bpsplit_by_row <- function(x, ncores, njobs=NULL) {
    if (ncores==1L) {
        return(list(x))
    }

    if (is.null(njobs)) {
        njobs <- bpnjobs_by_row(x, ncores)
    }

    last <- 0L
    out <- vector("list", ncores)
    for (i in seq_len(ncores)) {
        I <- last + seq_len(njobs[i])
        out[[i]] <- x[I,,drop=FALSE]
        last <- last + njobs[i]
    }

    out
}

bpsplit_by_col <- function(x, ncores, njobs=NULL) {
    if (ncores==1L) {
        return(list(x))
    }

    if (is.null(njobs)) {
        njobs <- bpnjobs_by_col(x, ncores)
    }

    last <- 0L
    out <- vector("list", ncores)
    for (i in seq_len(ncores)) {
        I <- last + seq_len(njobs[i])
        out[[i]] <- x[,I,drop=FALSE]
        last <- last + njobs[i]
    }

    out
}

.deviation <- function(njobs) {
    if (length(njobs)==1L) { return(0) }
    mean(abs(njobs - mean(njobs)))
}
