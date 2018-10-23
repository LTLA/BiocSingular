\name{runPCA}
\alias{runPCA}

\title{Principal components analysis}
\description{Perform a principal components analysis (PCA) on a target matrix with a specified SVD algorithm.}

\usage{
runPCA(x, rank, center=TRUE, scale=FALSE, BSPARAM=NULL, get.rotation=TRUE, 
    get.pcs=TRUE, ...) 
}

\arguments{
\item{x}{A numeric matrix-like object.}
\item{rank}{Integer scalar specifying the number of principal components to retain.}
\item{center}{A logical scalar indicating whether columns of \code{x} should be centered before the PCA is performed.
Alternatively, a numeric vector of length \code{ncol(x)} containing the value to subtract from each column of \code{x}.}
\item{scale}{A logical scalar indicating whether columns of \code{x} should be scaled to unit variance before the PCA is performed.
Alternatively, a numeric vector of length \code{ncol(x)} containing the scaling factor for each column of \code{x}.}
\item{BSPARAM}{A \linkS4class{BiocSingularParam} object specifying the algorithm to use.
This defaults to a \code{\link{ExactParam}} object when \code{NULL}.}
\item{get.rotation}{A logical scalar indicating whether rotation vectors should be returned.}
\item{get.pcs}{A logical scalar indicating whether the principal component scores should be returned.}
\item{...}{Further arguments to pass to \code{\link{runSVD}}, including \code{BPPARAM} to control parallelization.}
}

\details{
Any \code{scale} is applied after \code{center}, and mimics the behaviour of the \code{\link{scale}} function in \pkg{base}.
}

\value{
A list is returned containing:
\itemize{
\item{\code{sdev}, a numeric vector of length \code{rank} containing the standard deviations of the first \code{rank} principal components.}
\item{\code{rotation}, a numeric matrix with \code{rank} columns and \code{nrow(x)} rows, containing the first \code{rank} rotation vectors.
This is only returned if \code{get.rotation=TRUE}.}
\item{\code{x}, a numeric matrix with \code{rank} columns and \code{ncol(x)} rows, containing the scores for the first \code{rank} principal components.
This is only returned if \code{get.pcs=TRUE}.}
}
}

\author{
Aaron Lun
}

\seealso{
\code{\link{runSVD}} for the underlying SVD function.

\code{\link{BiocSingularParam}} for details on the algorithm choices.
}

\examples{
a <- matrix(rnorm(100000), ncol=20)
str(out <- runPCA(a, k=10))
}

