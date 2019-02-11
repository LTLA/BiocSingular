\name{runPCA}
\alias{runPCA}
\alias{runPCA,ANY-method}

\title{Principal components analysis}
\description{Perform a principal components analysis (PCA) on a target matrix with a specified SVD algorithm.}

\usage{
runPCA(x, ...)

\S4method{runPCA}{ANY}(x, rank, center=TRUE, scale=FALSE, get.rotation=TRUE,
    get.pcs=TRUE, ...) 
}

\arguments{
\item{x}{A numeric matrix-like object with samples as rows and variables as columns.}
\item{rank}{Integer scalar specifying the number of principal components to retain.}
\item{center}{A logical scalar indicating whether columns of \code{x} should be centered before the PCA is performed.
Alternatively, a numeric vector of length \code{ncol(x)} containing the value to subtract from each column of \code{x}.}
\item{scale}{A logical scalar indicating whether columns of \code{x} should be scaled to unit variance before the PCA is performed.
Alternatively, a numeric vector of length \code{ncol(x)} containing the scaling factor for each column of \code{x}.}
\item{get.rotation}{A logical scalar indicating whether rotation vectors should be returned.}
\item{get.pcs}{A logical scalar indicating whether the principal component scores should be returned.}
\item{...}{
For the generic, this contains arguments to pass to methods upon dispatch.

For the \code{ANY} method, this contains further arguments to pass to \code{\link{runSVD}}.
This includes \code{BSPARAM} to specify the algorithm that should be used, and \code{BPPARAM} to control parallelization.
}
}

\details{
This function simply calls \code{\link{runSVD}} and converts the results into a format similar to that returned by \code{\link{prcomp}}.

The generic is exported to allow other packages to implement their own \code{runPCA} methods for other \code{x} objects, e.g., \pkg{scater} for SingleCellExperiment inputs.
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

\code{?\linkS4class{BiocSingularParam}} for details on the algorithm choices.
}

\examples{
a <- matrix(rnorm(100000), ncol=20)
str(out <- runPCA(a, rank=10))
}

