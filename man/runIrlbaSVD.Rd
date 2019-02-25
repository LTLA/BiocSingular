\name{runIrlbaSVD}
\alias{runIrlbaSVD}

\title{Approximate SVD with \pkg{irlba}}
\description{Perform an approximate singular value decomposition with the augmented implicitly restarted Lanczos bidiagonalization algorithm.}

\usage{
runIrlbaSVD(x, k=5, nu=k, nv=k, center=FALSE, scale=FALSE, deferred=FALSE,
    extra.work=7, ..., fold=Inf, BPPARAM=SerialParam())
}

\arguments{
    \item{x}{A numeric matrix-like object to use in the SVD.}
    \item{k}{Integer scalar specifying the number of singular values to return.}
    \item{nu}{Integer scalar specifying the number of left singular vectors to return.}
    \item{nv}{Integer scalar specifying the number of right singular vectors to return.}
    \item{center}{A logical scalar indicating whether columns should be centered.
        Alternatively, a numeric vector or \code{NULL} - see \code{?"\link{BiocSingular-options}"}.} 
    \item{scale}{A logical scalar indicating whether columns should be scaled. 
        Alternatively, a numeric vector or \code{NULL} - see \code{?"\link{BiocSingular-options}"}.} 
    \item{deferred}{Logical scalar indicating whether centering/scaling should be deferred, see \code{?"\link{BiocSingular-options}"}.}
    \item{extra.work}{Integer scalar specifying the additional number of dimensions to use for the working subspace.}
    \item{...}{Further arguments to pass to \code{\link{irlba}}.}
    \item{fold}{Numeric scalar specifying the minimum fold difference between dimensions of \code{x} to compute the cross-product, 
        see \code{?"\link{BiocSingular-options}"}.}
    \item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.}
}

\details{
If \code{BPPARAM} has only 1 worker and a cross-product is not being computed, this function will use \code{\link{irlba}}'s own \code{center} and \code{scale} arguments.
This is effectively equivalent to deferred centering and scaling, despite the setting of \code{deferred=FALSE}.

For multiple workers, this function will parallelize all multiplication operations involving \code{x} according to the supplied \code{BPPARAM}.

The total dimensionality of the working subspace is defined as the maximum of \code{k}, \code{nu} and \code{nv}, plus the \code{extra.work}.
}

\value{
A list containing:
\itemize{
    \item \code{d}, a numeric vector of the first \code{k} singular values.
    \item \code{u}, a numeric matrix with \code{nrow(x)} rows and \code{nu} columns.
    Each column contains a left singular vector.
    \item \code{u}, a numeric matrix with \code{ncol(x)} rows and \code{nv} columns.
    Each column contains a right singular vector.
}
}

\author{
Aaron Lun
}

\seealso{
\code{\link{irlba}} for the underlying algorithm.
}

\examples{
a <- matrix(rnorm(100000), ncol=20)
out <- runIrlbaSVD(a)
str(out)
}
