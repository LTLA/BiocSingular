\name{runIrlba}
\alias{runIrlba}

\title{Approximate SVD with \pkg{irlba}}
\description{Perform an approximate singular value decomposition with \pkg{irlba}.}

\usage{
runIrlba(x, k=5, nu=k, nv=k, center=NULL, scale=NULL, extra.work=7, 
    ..., BPPARAM=SerialParam())
}

\arguments{
    \item{x}{A numeric matrix-like object to use in the SVD.}
    \item{k}{Integer scalar specifying the number of singular values to return.}
    \item{nu}{Integer scalar specifying the number of left singular vectors to return.}
    \item{nv}{Integer scalar specifying the number of right singular vectors to return.}
    \item{center}{Numeric vector of length equal to \code{ncol(x)}, specifying values to subtract from each column of \code{x}.}
    \item{scale}{Numeric vector of length equal to \code{ncol(x)}, specifying values to divide each column of \code{x}.
    This is applied after \code{center} if both are specified.}
    \item{extra.work}{Integer scalar specifying the additional number of dimensions to use for the working subspace.}
    \item{...}{Further arguments to pass to \code{\link{irlba}}.}
    \item{fold}{Numeric scalar specifying the minimum fold difference between dimensions of \code{x}, 
        at which the SVD is performed on the cross-product rather than \code{x} itself.}
    \item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.}
}

\details{
If \code{BPPARAM} has only 1 core \emph{and} \code{x} is an ordinary matrix or \code{dgCMatrix}, this function will use \code{\link{irlba}} with the default \code{fastpath=TRUE}.
Otherwise, it will use the slower R-based algorithm to take advantage of the generic (and parallelized) multiplication.

The total dimensionality of the working subspace is defined as the maximum of \code{k}, \code{nu} and \code{nv}, plus the \code{extra.work}.

If the dimensions of \code{x} are very different (as defined by \code{fold}), \code{\link{irlba}} is called on the cross-product of \code{x} instead.
This can be much faster than directly performing the SVD on a very fat or tall \code{x}.
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
out <- runIrlba(a)
str(out)
}
