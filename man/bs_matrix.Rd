\name{Internal methods}
\alias{deferred-operations}
\alias{\%*\%,ANY,bs_matrix-method}
\alias{\%*\%,bs_matrix,ANY-method}
\alias{\%*\%,bs_matrix,bs_matrix-method}
\alias{dim,bs_matrix-method}
\alias{dimnames,bs_matrix-method}
\alias{length,bs_matrix-method}
\docType{class}

\title{Deferred centering/scaling}
\description{An overview of how centering and scaling can be deferred during matrix multiplication and cross-product calculations.}

\details{
The \pkg{BiocSingular} package has a specialized matrix class that defers centering (and to some extent, scaling) during matrix multiplication.
The matrix multiplication is performed on the original matrix, and then the centering/scaling operations are applied to the matrix product.
This allows direct use of the \code{\%*\%} method for each matrix representation, to exploit features of the underlying matrix representation for greater speed.

Unfortunately, the speed-up with deferred centering comes at the cost of increasing the risk of catastrophic cancellation.
The procedure requires subtraction of one large intermediate number from another to obtain the values of the final matrix product.
This could result in a loss of numerical precision that compromises the accuracy of the various SVD algorithms. 

The default approach is to explicitly create a dense in-memory centred/scaled matrix via block processing (see \code{\link{blockGrid}} in the \pkg{DelayedArray} package).
This avoids problems with numerical precision as large intermediate values are not formed.
In doing so, we consistently favour accuracy over speed unless the functions are specifically instructed to do otherwise, i.e., with \code{deferred=TRUE}.
}

\author{
Aaron Lun
}
