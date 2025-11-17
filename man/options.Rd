\name{BiocSingular options}
\alias{BiocSingular-options}

\title{Global SVD options}
\description{An overview of the available options when performing SVD with any algorithm.}

\section{Computing the cross-product}{ 
If the dimensions of the input matrix are very different, it may be faster to compute the cross-product and perform the SVD on the resulting square matrix,
rather than performing SVD directly on a very fat or tall input matrix.
The cross-product can often be computed very quickly due to good data locality, yielding a small square matrix that is easily handled by any SVD algorithm.
This is especially true in cases where the input matrix is not held in memory.
Calculation of the cross-product only involves one read across the entire data set, 
while direct application of approximate methods like \code{\link[irlba]{irlba}} or \code{\link[rsvd]{rsvd}} would need to access the data multiple times.

The various \pkg{BiocSingular} SVD functions allow users to specify the minimum fold difference (via the \code{fold} argument) at which a cross-product should be computed. 
Setting \code{fold=1} will always compute the cross-product for any matrix - this is probably unwise. 
By contrast, setting \code{fold=Inf} means that the cross-product is never computed.
This is currently the default in all functions, to provide the most expected behaviour unless specifically instructed otherwise.
}

\section{Centering and scaling}{
In general, each SVD function performs the SVD on \code{t((t(x) - C)/S)} where \code{C} and \code{S} are numeric vectors of length equal to \code{ncol(x)}.
The values of \code{C} and \code{S} are defined according to the \code{center} and \code{scale} options.
\itemize{
    \item If \code{center=TRUE}, \code{C} is defined as the column means of \code{x}.
        If \code{center=NULL} or \code{FALSE}, all elements of \code{C} are set to zero.
        If \code{center} is a numeric vector with length equal to \code{ncol(x)}, it is used to directly define \code{C}.
    \item If \code{scale=TRUE}, the \code{i}th element of \code{S} is defined as the square root of \code{sum((x[,i] - C[i])^2)/(ncol(x)-1)}, for whatever \code{C} was defined above.
        This mimics the behaviour of \code{\link{scale}}.
        If \code{scale=NULL} or \code{FALSE}, all elements of \code{S} are set to unity.
        If \code{scale} is a numeric vector with length equal to \code{ncol(x)}, it is used to directly define \code{S}.
}

Setting \code{center} or \code{scale} is more memory-efficient than modifiying the input \code{x} directly.
This is because the function will avoid constructing intermediate centered (possibly non-sparse) matrices.
}

\section{Deferred centering and scaling}{
Many of the SVD algorithms (and computation of the cross-product) involve repeated matrix multiplications.
We speed this up by using the \link[ScaledMatrix]{ScaledMatrix} class to defer centering (and to some extent, scaling) during matrix multiplication.
The matrix multiplication is performed on the original matrix, and then the centering/scaling operations are applied to the matrix product.
This allows direct use of the \code{\%*\%} method for each matrix representation, to exploit features of the underlying representation (e.g., sparsity) for greater speed.

Unfortunately, the speed-up with deferred centering comes at the cost of increasing the risk of catastrophic cancellation.
The procedure requires subtraction of one large intermediate number from another to obtain the values of the final matrix product.
This could result in a loss of numerical precision that compromises the accuracy of the various SVD algorithms. 

The default approach is to explicitly create a dense in-memory centred/scaled matrix,
possibly via block processing (see \code{\link[DelayedArray]{defaultAutoGrid}} from the \pkg{DelayedArray} package).
This avoids problems with numerical precision as large intermediate values are not formed.
In doing so, we consistently favour accuracy over speed unless the functions are specifically instructed to do otherwise, i.e., with \code{deferred=TRUE}.
}

\author{
Aaron Lun
}
