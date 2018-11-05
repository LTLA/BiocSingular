\name{Internal methods}
\alias{\%*\%,ANY,bs_matrix-method}
\alias{\%*\%,bs_matrix,ANY-method}
\alias{\%*\%,bs_matrix,bs_matrix-method}
\alias{dim,bs_matrix-method}
\alias{dimnames,bs_matrix-method}
\alias{length,bs_matrix-method}
\docType{class}

\title{Internal S4 methods}
\description{Internal S4 methods for the \pkg{BiocSingular} matrix class.
This should not be relevant to most users.} 

\details{
The \pkg{BiocSingular} package has a specialized matrix class that defers centering (and to some extent, scaling) during matrix multiplication.
This allows multiplication methods to exploit features of the underlying matrix representation for greater speed.
Otherwise, explicit coercion to a dense in-memory matrix would be required prior to multiplication.

The cost of greater speed is the risk of catastrophic cancellation when large numbers are added and subtracted.
This could result in a loss of numerical precision that compromises the accuracy of the various SVD algorithms. 
In practice, this does not seem to be a major concern for algorithms using double-precision values.

One can check the effect on precision by using a \linkS4class{DelayedArray} with delayed centering and scaling operations. 
This will be explicitly realized into a dense array during matrix multiplication, avoiding the need for subtraction with deferred centering.
The result should be more numerically stable and can be compared to the result of using \code{center} and \code{scale}.
}

\author{
Aaron Lun
}
