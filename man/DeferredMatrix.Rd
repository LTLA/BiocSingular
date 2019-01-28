\name{DeferredMatrix}

\alias{DeferredMatrix}
\alias{DeferredMatrix-class}

\alias{dim,DeferredMatrix-method}
\alias{dimnames,DeferredMatrix-method}
\alias{dimnames<-,DeferredMatrix,ANY-method}
\alias{length,DeferredMatrix-method}
\alias{show,DeferredMatrix-method}

\alias{as.matrix,DeferredMatrix-method}
\alias{t,DeferredMatrix-method}

\alias{[,DeferredMatrix,ANY,ANY,ANY-method}

\alias{colSums,DeferredMatrix-method}
\alias{rowSums,DeferredMatrix-method}

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

\alias{extract_array,DeferredMatrix-method}

\docType{class}

\title{The DeferredMatrix class}
\description{
Methods for the DeferredMatrix class, supporting deferred centering and scaling.
}

\usage{
DeferredMatrix(x, center=NULL, scale=NULL)
}

\arguments{
\item{x}{A matrix-like object.}
\item{center}{A numeric vector of length equal to \code{ncol(x)}, where each element is to be subtracted from the corresponding column of \code{x}.
A \code{NULL} value indicates that no subtraction is to be performed.}
\item{scale}{A numeric vector of length equal to \code{ncol(x)}, where each element is to divide from the corresponding column of \code{x} (after subtraction).
A \code{NULL} value indicates that no division is to be performed.}
}

\value{
A DeferredMatrix object equivalent to \code{t((t(x) - center)/scale)}.
}

\section{Available methods}{
DeferredMatrix objects support standard operations like \code{dim}, \code{dimnames}, \code{dimnames<-}, \code{length} and \code{show}.

Calling \code{as.matrix} will return an ordinary R matrix with centering and scaling operations applied.

Calling \code{t} will indicate that the matrix is to be transposed, along with any scaling/centering operations.
This can effectively be used to allow centering/scaling on the rows if the input \code{x} is transposed.

Matrix subsetting will return a DeferredMatrix of the specified dimensions.

Matrix multiplication via \code{\%*\%}, \code{crossprod} and \code{tcrossprod} are supported.
Multiplying two DeferredMatrix objects is not yet supported.

Calculation of row and column sums is supported by \code{colSums} and \code{rowSums}.
Other matrix statistics are not yet supported.

DeferredMatrix objects can be wrapped in \code{\link{DelayedArray}} instances to perform otherwise unsupported calculations via block processing.
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
