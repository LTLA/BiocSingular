\name{runExactSVD}
\alias{runExactSVD}

\title{Exact SVD}
\description{Perform an exact singular value decomposition.}

\usage{
runExactSVD(x, k=min(dim(x)), nu=k, nv=k, center=NULL, scale=NULL, fold=5)
}

\arguments{
    \item{x}{A numeric matrix-like object to use in the SVD.}
    \item{k}{Integer scalar specifying the number of singular values to return.}
    \item{nu}{Integer scalar specifying the number of left singular vectors to return.}
    \item{nv}{Integer scalar specifying the number of right singular vectors to return.}
    \item{center}{Numeric vector of length equal to \code{ncol(x)}, specifying values to subtract from each column of \code{x}.}
    \item{scale}{Numeric vector of length equal to \code{ncol(x)}, specifying values to divide each column of \code{x}.
    This is applied after \code{center} if both are specified.}
    \item{fold}{Numeric scalar specifying the minimum fold difference between dimensions of \code{x}, 
        at which the SVD is performed on the cross-product rather than \code{x} itself.}
}
        
\details{
If the dimensions of \code{x} are very different, a cross-product of \code{x} is computed and the SVD is performed on the resulting square matrix.
This can be much faster than directly performing the SVD on a very fat or tall \code{x}.

Supplying \code{center} or \code{scale} may be more efficient than directly scaling \code{x} when computing a cross-product.
This is because the function will avoid constructing intermediate centered (possibly non-sparse) matrices.
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