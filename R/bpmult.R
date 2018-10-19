#' @importFrom BiocParallel bplapply bpnworkers
bpmult <- function(x, y, BPPARAM) 
# Parallelizing the multiplication by splitting the work by rows of 'x' or columns of 'y',
# depending on which one yields a smaller per-worker serialization.
{ 
    ncores <- bpnworkers(BPPARAM)
    y.not.mat <- is.null(dim(y))
    x.not.mat <- is.null(dim(x))
    if (ncores==1L || (x.not.mat && y.not.mat)) {
        return(x %*% y)

    } else if (y.not.mat || object.size(x) > object.size(y)) { # 'x' is larger, so we break it up.
        
        last <- 0L
        inc <- ceiling(nrow(x)/ncores)
        left <- vector("list", ncores)

        for (i in seq_len(ncores)) {
            I <- last + seq_len(inc)
            I <- I[I <= nrow(x)]
            left[[i]] <- x[I,,drop=FALSE]
            last <- last + inc
        }

        out <- bplapply(FUN=.internal_mult, left, right=y, BPPARAM=BPPARAM)
        return(do.call(rbind, out))

    } else { # 'y' is larger, so we break it up instead.
        last <- 0L
        inc <- ceiling(ncol(y)/ncores)
        right <- vector("list", ncores)

        for (i in seq_len(ncores)) {
            I <- last + seq_len(inc)
            I <- I[I <= ncol(y)]
            right[[i]] <- y[,I,drop=FALSE]
            last <- last + inc
        }

        out <- bplapply(FUN=.internal_mult, right, left=x, BPPARAM=BPPARAM)
        return(do.call(cbind, out))

    }
}

.internal_mult <- function(left, right) 
# Avoid problems with using 'y' in lapply on '%*%' directly.
{
    left %*% right
}
