\name{DeferredMatrix}

\alias{DeferredMatrixSeed}
\alias{DeferredMatrixSeed-class}

\alias{dim,DeferredMatrixSeed-method}
\alias{dimnames,DeferredMatrixSeed-method}
\alias{extract_array,DeferredMatrixSeed-method}
\alias{DelayedArray,DeferredMatrixSeed-method}
\alias{show,DeferredMatrixSeed-method}

\alias{DeferredMatrix}
\alias{DeferredMatrix-class}

\alias{dimnames<-,DeferredMatrix,ANY-method}
\alias{t,DeferredMatrix-method}
\alias{[,DeferredMatrix,ANY,ANY,ANY-method}

\alias{colSums,DeferredMatrix-method}
\alias{rowSums,DeferredMatrix-method}
\alias{colMeans,DeferredMatrix-method}
\alias{rowMeans,DeferredMatrix-method}

\alias{\%*\%,ANY,DeferredMatrix-method}
\alias{\%*\%,DeferredMatrix,ANY-method}
\alias{\%*\%,DeferredMatrix,DeferredMatrix-method}

\alias{crossprod,DeferredMatrix,missing-method}
\alias{crossprod,DeferredMatrix,ANY-method}
\alias{crossprod,ANY,DeferredMatrix-method}
\alias{crossprod,DeferredMatrix,DeferredMatrix-method}

\alias{tcrossprod,DeferredMatrix,missing-method}
\alias{tcrossprod,DeferredMatrix,ANY-method}
\alias{tcrossprod,ANY,DeferredMatrix-method}
\alias{tcrossprod,DeferredMatrix,DeferredMatrix-method}

\docType{class}

\title{The DeferredMatrix class}
\description{
Definitions of the DeferredMatrixSeed and DeferredMatrix classes and their associated methods.
These classes are designed to support deferred centering and scaling of the columns prior to a principal components analysis.
}

\usage{
DeferredMatrixSeed(x, center=NULL, scale=NULL)

DeferredMatrix(x, center=NULL, scale=NULL)
}

\arguments{
\item{x}{A matrix-like object.

This can alternatively be a DeferredMatrixSeed, in which case any values of \code{center} and \code{scale} are ignored.}
\item{center}{A numeric vector of length equal to \code{ncol(x)}, where each element is to be subtracted from the corresponding column of \code{x}.
A \code{NULL} value indicates that no subtraction is to be performed.}
\item{scale}{A numeric vector of length equal to \code{ncol(x)}, where each element is to divide from the corresponding column of \code{x} (after subtraction).
A \code{NULL} value indicates that no division is to be performed.}
}

\value{
The \code{DeferredMatrixSeed} constructor will return a DeferredMatrixSeed object.

The \code{DeferredMatrix} constructor will return a DeferredMatrix object equivalent to \code{t((t(x) - center)/scale)}.
}

\section{Methods for DeferredMatrixSeed objects}{
DeferredMatrixSeed objects are implemented as \linkS4class{DelayedMatrix} backends.
They support standard operations like \code{dim}, \code{dimnames} and \code{extract_array}.

Passing a DeferredMatrixSeed object to the \code{\link{DelayedArray}} constructor will create a DeferredMatrix object.

It is possible for \code{x} to contain a DeferredMatrix, thus nesting one DeferredMatrix inside another.
This can occasionally be useful in combination with transposition to achieve centering/scaling in both dimensions.
}

\section{Methods for DeferredMatrix objects}{
DeferredMatrix objects are derived from \linkS4class{DelayedMatrix} objects and support all of valid operations on the latter.
Several functions are specialized for greater efficiency when operating on DeferredMatrix instances, including:
\itemize{
    \item Subsetting, transposition and replacement of row/column names.
        These will return a new DeferredMatrix rather than a DelayedMatrix.
    \item Matrix multiplication via \code{\%*\%}, \code{crossprod} and \code{tcrossprod}.
        These functions will return a DelayedMatrix.
        Also see \code{?"\link{BiocSingular-options}"}.
    \item Calculation of row and column sums and means by \code{colSums}, \code{rowSums}, etc. 
}

All other operations applied to a DeferredMatrix will use the underlying \pkg{DelayedArray} machinery.
Unary or binary operations will generally create a new DelayedMatrix instance containing a DeferredMatrixSeed.

Tranposition can effectively be used to allow centering/scaling on the rows if the input \code{x} is transposed.
}

\author{
Aaron Lun
}

\examples{
library(Matrix)
y <- DeferredMatrix(rsparsematrix(10, 20, 0.1), 
    center=rnorm(20), scale=1+runif(20))
y

crossprod(y)
tcrossprod(y)
y \%*\% rnorm(20)
}

\seealso{
\code{?"\link{BiocSingular-options}"} for comments about numerical precision with deferred centering and scaling.
}
