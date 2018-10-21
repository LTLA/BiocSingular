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

bpsplit_by_row <- function(x, njobs) {
    last <- 0L
    out <- vector("list", length(njobs))
    for (i in seq_along(njobs)) {
        I <- last + seq_len(njobs[i])
        out[[i]] <- x[I,,drop=FALSE]
        last <- last + njobs[i]
    }
    out
}

bpsplit_by_col <- function(x, njobs) {
    last <- 0L
    out <- vector("list", length(njobs))
    for (i in seq_along(njobs)) {
        I <- last + seq_len(njobs[i])
        out[[i]] <- x[,I,drop=FALSE]
        last <- last + njobs[i]
    }
    out
}

dispatcher <- function(x, y, ncores, x.transposed=FALSE, y.transposed=FALSE) 
# Chooses how to dispatch jobs for matrix-multiplication-like operations.
# The aim is to minimize the maximum amount of data that needs to be sent
# to any single core. This reduces I/O operations per core and coincides
# with reducing the maximum amount of calculations for each core.
{
    x.by.row <- bpnjobs_by_row(x, ncores)
    x.by.col <- bpnjobs_by_col(x, ncores)
    y.by.row <- bpnjobs_by_row(y, ncores)
    y.by.col <- bpnjobs_by_col(y, ncores)

    # Determining the parallelization choices.
    if (!x.transposed && !y.transposed) { # %*%
        x.choice <- x.by.row
        y.choice <- y.by.col
        common.dim <- ncol(x)

        both.choice.x <- x.by.col
        x.other.dim <- nrow(x)
        both.choice.y <- y.by.row
        y.other.dim <- ncol(y)

    } else if (x.transposed && !y.transposed) { # crossprod
        x.choice <- x.by.col
        y.choice <- y.by.col
        common.dim <- nrow(x)

        both.choice.x <- x.by.row
        x.other.dim <- ncol(x)
        both.choice.y <- y.by.row
        y.other.dim <- ncol(y)

    } else if (y.transposed && !x.transposed) { # tcrossprod
        x.choice <- x.by.row
        y.choice <- y.by.row
        common.dim <- ncol(x)

        both.choice.x <- x.by.col
        x.other.dim <- nrow(x)
        both.choice.y <- y.by.col
        y.other.dim <- nrow(y)

    } else {
        stop("not supported")
    }

    # Computing the amount of data sent under each scheme.
    # This is not quite accurate for alternative representations or chunking layouts,
    # nor for various parallelization schemes (depending on whether serialization is necessary).
    # However, refining the calculations increases the complexity considerably,
    # and it should be pretty clear which one is better when one matrix is large and file-backed.
    x.cost <- x.choice * common.dim + nrow(y) * ncol(y)
    y.cost <- y.choice * common.dim + nrow(x) * ncol(x)
    both.cost.x <- both.choice.x * (x.other.dim + y.other.dim)
    both.cost.y <- both.choice.y * (x.other.dim + y.other.dim)

    collected <- c(x=max(x.cost), y=max(y.cost), both.x=max(both.cost.x), both.y=max(both.cost.y))
    best <- names(collected)[which.min(collected)]
    switch(best,
        x=list(choice="x", jobs=x.choice),
        y=list(choice="y", jobs=y.choice),
        both.x=list(choice="both", jobs=both.choice.x),
        both.y=list(choice="both", jobs=both.choice.y)
    )
}
