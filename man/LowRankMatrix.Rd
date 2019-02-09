\name{LowRankMatrix}

\alias{LowRankMatrixSeed}
\alias{LowRankMatrixSeed-class}

\alias{dim,LowRankMatrixSeed-method}
\alias{dimnames,LowRankMatrixSeed-method}
\alias{extract_array,LowRankMatrixSeed-method}
\alias{DelayedArray,LowRankMatrixSeed-method}
\alias{show,LowRankMatrixSeed-method}

\alias{LowRankMatrix}
\alias{LowRankMatrix-class}

\alias{dimnames<-,LowRankMatrix,ANY-method}
\alias{t,LowRankMatrix-method}
\alias{[,LowRankMatrix,ANY,ANY,ANY-method}

\docType{class}

\title{The LowRankMatrix class}
\description{
Definitions of the LowRankMatrixSeed and LowRankMatrix classes and their associated methods.
These classes are designed to provide a memory-efficient representation of a low-rank reconstruction, e.g., after a principal components analysis.
}

\usage{
LowRankMatrixSeed(rotation, components)

LowRankMatrix(rotation, components)
}

\arguments{
\item{rotation}{A matrix-like object where each row corresponds to a row of the LowRankMatrix object.

For the \code{LowRankMatrix} constructor, this can alternatively be a LowRankMatrixSeed, in which case \code{components} does not need to be specified.}
\item{components}{A matrix-like object where each row corresponds to a column of the LowRankMatrix object.}
}

\value{
The \code{LowRankMatrixSeed} constructor will return a LowRankMatrixSeed object.

The \code{LowRankMatrix} constructor will return a LowRankMatrix object equivalent to \code{tcrossprod(rotation, components)}.
}

\section{Methods for LowRankMatrixSeed objects}{
LowRankMatrixSeed objects are implemented as \linkS4class{DelayedMatrix} backends.
They support standard operations like \code{dim}, \code{dimnames} and \code{extract_array}.

Passing a LowRankMatrixSeed object to the \code{\link{DelayedArray}} constructor will create a LowRankMatrix object.
}

\section{Methods for LowRankMatrix objects}{
LowRankMatrix objects are derived from \linkS4class{DelayedMatrix} objects and support all of valid operations on the latter.
Subsetting, transposition and replacement of row/column names are specialized for greater efficiency when operating on LowRankMatrix instances,
and will return a new LowRankMatrix rather than a DelayedMatrix.

All other operations applied to a LowRankMatrix will use the underlying \pkg{DelayedArray} machinery.
Unary or binary operations will generally create a new DelayedMatrix instance containing a LowRankMatrixSeed.
}

\author{
Aaron Lun
}

\examples{
a <- matrix(rnorm(100000), ncol=20)
out <- runPCA(a, rank=10)

lr <- LowRankMatrix(out$rotation, out$x)
}

\seealso{
\code{\link{runPCA}} to generate the rotation and component matrices.
}
