\name{runRandomSVD}
\alias{runRandomSVD}

\title{Approximate SVD with \pkg{rsvd}}
\description{Perform a randomized singular value decomposition.}

\usage{
runRandomSVD(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, deferred=FALSE,
    ..., fold=5, BPPARAM=SerialParam())
}

\arguments{
    \item{x}{A numeric matrix-like object to use in the SVD.}
    \item{k}{Integer scalar specifying the number of singular values to return.}
    \item{nu}{Integer scalar specifying the number of left singular vectors to return.}
    \item{nv}{Integer scalar specifying the number of right singular vectors to return.}
    \item{center}{Numeric vector of length equal to \code{ncol(x)}, specifying values to subtract from each column of \code{x}.}
    \item{scale}{Numeric vector of length equal to \code{ncol(x)}, specifying values to divide each column of \code{x}.
    This is applied after \code{center} if both are specified.}
    \item{deferred}{Logical scalar indicating whether centering/scaling should be deferred, see \code{?"\link{BiocSingular-options}"}.}
    \item{...}{Further arguments to pass to \code{\link{rsvd}}.}
    \item{fold}{Numeric scalar specifying the minimum fold difference between dimensions of \code{x} to compute the cross-product, 
        see \code{?"\link{BiocSingular-options}"}.}
    \item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.}
}

\details{
Supplying \code{center} or \code{scale} is more efficient than directly scaling \code{x} when computing a cross-product.
This is because the function will avoid constructing intermediate centered (possibly non-sparse) matrices.
        
All multiplication operations in \code{\link{rsvd}} involving \code{x} will be parallelized according to the supplied \code{BPPARAM}.

The dimensionality of the working subspace is defined as the maximum of \code{k}, \code{nu} and \code{nv}, plus the \code{q} specified in \code{...}.
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
\code{\link{rsvd}} for the underlying algorithm.
}

\examples{
a <- matrix(rnorm(100000), ncol=20)
out <- runRandomSVD(a)
str(out)
}
