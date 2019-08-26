\name{ResidualMatrix}

\alias{ResidualMatrixSeed}
\alias{ResidualMatrixSeed-class}

\alias{dim,ResidualMatrixSeed-method}
\alias{dimnames,ResidualMatrixSeed-method}
\alias{extract_array,ResidualMatrixSeed-method}
\alias{DelayedArray,ResidualMatrixSeed-method}
\alias{show,ResidualMatrixSeed-method}

\alias{ResidualMatrix}
\alias{ResidualMatrix-class}

\alias{dimnames<-,ResidualMatrix,ANY-method}
\alias{t,ResidualMatrix-method}
\alias{[,ResidualMatrix,ANY,ANY,ANY-method}

\alias{colSums,ResidualMatrix-method}
\alias{rowSums,ResidualMatrix-method}
\alias{colMeans,ResidualMatrix-method}
\alias{rowMeans,ResidualMatrix-method}

\alias{\%*\%,ANY,ResidualMatrix-method}
\alias{\%*\%,ResidualMatrix,ANY-method}
\alias{\%*\%,ResidualMatrix,ResidualMatrix-method}

\alias{crossprod,ResidualMatrix,missing-method}
\alias{crossprod,ResidualMatrix,ANY-method}
\alias{crossprod,ANY,ResidualMatrix-method}
\alias{crossprod,ResidualMatrix,ResidualMatrix-method}

\alias{tcrossprod,ResidualMatrix,missing-method}
\alias{tcrossprod,ResidualMatrix,ANY-method}
\alias{tcrossprod,ANY,ResidualMatrix-method}
\alias{tcrossprod,ResidualMatrix,ResidualMatrix-method}

\alias{runPCA,ResidualMatrix-method}

\docType{class}

\title{The ResidualMatrix class}
\description{
Definitions of the ResidualMatrixSeed and ResidualMatrix classes and their associated methods.
These classes are designed to support delayed calculation of the residuals from a linear model fit, 
usually prior to a principal components analysis.
The aim is to perform matrix multiplication without explicitly calculating the residuals,
allowing efficient computation based on features of the the original matrix (e.g., sparsity).
}

\usage{
ResidualMatrixSeed(x, design=NULL)

ResidualMatrix(x, design=NULL)
}

\arguments{
\item{x}{A matrix-like object.

This can alternatively be a ResidualMatrixSeed, in which case \code{design} is ignored.}
\item{design}{A numeric matrix containing the experimental design,
to be used for linear model fitting on each \emph{column} of \code{x}.}
}

\value{
The \code{ResidualMatrixSeed} constructor will return a ResidualMatrixSeed object.

The \code{ResidualMatrix} constructor will return a ResidualMatrix object,
containing values equivalent to \code{lm.fit(x=design, y=x)$residuals}.
}

\section{Methods for ResidualMatrixSeed objects}{
ResidualMatrixSeed objects are implemented as \linkS4class{DelayedMatrix} backends.
They support standard operations like \code{dim}, \code{dimnames} and \code{extract_array}.

Passing a ResidualMatrixSeed object to the \code{\link{DelayedArray}} constructor will create a ResidualMatrix object.
}

\section{Methods for ResidualMatrix objects}{
ResidualMatrix objects are derived from \linkS4class{DelayedMatrix} objects and support all of valid operations on the latter.
Several functions are specialized for greater efficiency when operating on ResidualMatrix instances, including:
\itemize{
    \item Subsetting, transposition and replacement of row/column names.
        These will return a new ResidualMatrix rather than a DelayedMatrix.
    \item Matrix multiplication via \code{\%*\%}, \code{crossprod} and \code{tcrossprod}.
        These functions will return a DelayedMatrix.
    \item Calculation of row and column sums and means by \code{colSums}, \code{rowSums}, etc. 
}

All other operations applied to a ResidualMatrix will use the underlying \pkg{DelayedArray} machinery.
Unary or binary operations will generally create a new DelayedMatrix instance containing a ResidualMatrixSeed.
}

\section{PCA with ResidualMatrix objects}{
\code{runPCA(x, center=TRUE, ...)}
will perform a PCA on a ResidualMatrix object \code{x},
with arguments in \code{...} passed to the ANY method for \code{\link{runPCA}}.
The \linkS4class{ResidualMatrix} method has the special behaviour that \code{center=TRUE} is ignored if:
\itemize{
    \item \code{x} was generated using a \code{design} that can be parameterized with an intercept.
    This means that the residuals on each column of \code{x} are centered at zero already.
    \item No subsetting by row was performed on \code{x},
    i.e., the zero-centering of the residuals is preserved.
}
This improves efficiency by avoiding an unnecessary additional centering step,
which would otherwise require block processing or deferred centering - see \code{?"\link{BiocSingular-options}"}.

The ResidualMatrix is particularly efficient when combined with approximate PCA strategies based on matrix multiplication.
This is achieved by setting \code{BSPARAM} to values like \code{\link{IrlbaParam}} or \code{\link{RandomParam}}.
The matrix product used in each algorithms can be computed efficiently without actually computing the residuals.
}

\author{
Aaron Lun
}

\examples{
design <- model.matrix(~gl(5, 50))

library(Matrix)
y0 <- rsparsematrix(nrow(design), 200, 0.1)
y <- ResidualMatrix(y0, design)
y

# For comparison.
fit <- lm.fit(x=design, y=as.matrix(y0))
DelayedArray(fit$residuals)

crossprod(y)
tcrossprod(y)
y \%*\% rnorm(200)

# PCA can be performed very quickly on ResidualMatrix 
# instances as the underlying representation can be 
# used directly, e.g., without loss of sparsity.
pc.out <- runPCA(y, 10, BSPARAM=IrlbaParam())
}
