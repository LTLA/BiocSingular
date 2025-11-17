\name{runExactSVD}
\alias{runExactSVD}

\title{Exact SVD}
\description{Perform an exact singular value decomposition.}

\usage{
runExactSVD(x, k=min(dim(x)), nu=k, nv=k, center=FALSE, scale=FALSE, 
    deferred=FALSE, fold=Inf, BPPARAM=SerialParam())
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
    \item{fold}{Numeric scalar specifying the minimum fold difference between dimensions of \code{x} to compute the cross-product, 
        see \code{?"\link{BiocSingular-options}"}.}
    \item{BPPARAM}{A \link[BiocParallel]{BiocParallelParam} object specifying how parallelization should be performed.}
}

\details{
If any of \code{k}, \code{nu} or \code{nv} exceeds \code{min(dim(x))}, they will be capped and a warning will be raised.
The exception is when they are explicitly set to \code{Inf}, in which case all singular values/vectors of \code{x} are returned without any warning.

Note that parallelization via \code{BPPARAM} is only applied to the calculation of the cross-product.
It has no effect for near-square matrices where the SVD is computed directly.
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
\code{\link{svd}} for the underlying algorithm.
}

\examples{
a <- matrix(rnorm(100000), ncol=20)
out <- runExactSVD(a)
str(out)
}
