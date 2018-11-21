\name{BiocSingular options}
\alias{BiocSingular-options}
\alias{\%*\%,ANY,bs_matrix-method}
\alias{\%*\%,bs_matrix,ANY-method}
\alias{\%*\%,bs_matrix,bs_matrix-method}
\alias{dim,bs_matrix-method}
\alias{dimnames,bs_matrix-method}
\alias{length,bs_matrix-method}
\docType{class}

\title{Global SVD options}
\description{An overview of the available options when performing SVD with any algorithm.}

\section{Computing the cross-product}{ 
If the dimensions of the input matrix are very different, it may be faster to compute the cross-product and perform the SVD on the resulting square matrix,
rather than performing SVD directly on a very fat or tall input matrix.
The cross-product can often be computed very quickly due to good data locality, yielding a small square matrix that is easily handled by any SVD algorithm.
This is especially true in cases where the input matrix is not held in memory.
Calculation of the cross-product only involves one read across the entire data set, 
while direct application of approximate methods like \code{\link{irlba}} or \code{\link{rsvd}} would need to access the data multiple times.

The various \pkg{BiocSingular} SVD functions allow users to specify the minimum fold difference (via the \code{fold} argument) at which a cross-product should be computed. 
The default is to require a difference of 5-fold between the two dimensions.
This avoids an unnecessary calculation of the cross-product when the matrix is already close to square-shaped.
Setting \code{fold=Inf} means that the cross-product is never computed, which may be useful when even the cross-product is too large to handle in memory.
Setting \code{fold=1} will always compute the cross-product - this is probably unwise. 
}

\section{Deferred centering and scaling}{
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
