\name{BiocSingularParam}
\alias{BiocSingularParam-class}

\alias{ExactParam}
\alias{ExactParam-class}

\alias{IrlbaParam}
\alias{IrlbaParam-class}
\alias{show,IrlbaParam-method}

\alias{RandomParam}
\alias{RandomParam-class}
\alias{show,RandomParam-method}

\alias{FastAutoParam}
\alias{FastAutoParam-class}

\docType{class}
\alias{show,BiocSingularParam-method}

\alias{bsfold}
\alias{bsdeferred}
\alias{bsparam}

\title{BiocSingularParam classes}
\description{Classes for specifying the type of singular value decomposition (SVD) algorithm and associated parameters.}

\usage{
ExactParam(deferred=FALSE, fold=Inf)

IrlbaParam(deferred=FALSE, fold=Inf, extra.work=7, ...)

RandomParam(deferred=FALSE, fold=Inf, ...)

FastAutoParam(deferred=FALSE, fold=Inf)

bsparam()
}

\arguments{
    \item{deferred}{Logical scalar indicating whether centering/scaling should be deferred, see \code{?"\link{BiocSingular-options}"}.}
    \item{fold}{Numeric scalar specifying the minimum fold-difference for cross-product calculation, see \code{?"\link{BiocSingular-options}"}.}
    \item{extra.work}{Integer scalar, additional dimensionality of the workspace in \code{\link{runIrlbaSVD}}.}
    \item{...}{Additional arguments to pass to \code{\link{runIrlbaSVD}} or \code{\link{runRandomSVD}}.
    This should not include any named arguments in those functions.}
}

\details{
The BiocSingularParam class controls dispatch of the \code{\link{runSVD}} generic to specific algorithms for SVD.
The BiocSingularParam class itself is virtual, with several concrete subclasses available in this package:
\describe{
    \item{}{\code{ExactParam}: exact SVD with \code{\link{runExactSVD}}.}
    \item{}{\code{IrlbaParam}: approximate SVD with \pkg{irlba} via \code{\link{runIrlbaSVD}}.}
    \item{}{\code{RandomParam}: approximate SVD with \pkg{rsvd} via \code{\link{runRandomSVD}}.}
    \item{}{\code{FastAutoParam}: fast approximate SVD, chosen based on the matrix representation.}
}
These objects also hold parameters specifying how each algorithm should be run on an arbitrary data set.
See the associated documentation pages for each specific SVD method for more details.
}

\section{Methods}{
In the following code snippets, \code{x} is a BiocSingularParam object.
\describe{
\item{\code{show(object)}:}{Display the class of \code{object} and a summary of the set parameters.}
\item{\code{bsfold(object)}:}{Return a numeric scalar specifying the fold-difference for cross-product calculation,
see \dQuote{Computing the cross-product} in \code{?"\link{BiocSingular-options}"}.}
\item{\code{bsdeferred(object)}:}{Return a logical scalar indicating whether centering and scaling should be deferred.
see \dQuote{Deferred centering and scaling} in \code{?"\link{BiocSingular-options}"}.}
}
}

\section{Automatically choosing a fast SVD algorithm}{
Running \code{\link{runSVD}} with a FastAutoParam object will dispatch to \code{\link{runIrlbaSVD}} for most matrices.
This is the default choice as IRLBA is fast and the approximation is highly similar to the exact SVD result. 
The exception is that of instances of the \linkS4class{DelayedMatrix} class, or any subclass that lacks its own specialized \code{\%*\%} method.
In such cases, \code{\link{runSVD}} with a FastAutoParam will dispatch to \code{\link{runRandomSVD}}, which minimizes the number of matrix multiplication steps and the associated costly block processing.
However, if a DelayedMatrix subclass has its own \code{\%*\%} method, it is assumed to be efficient enough to afford repeated multiplications in \code{\link{runIrlbaSVD}}.
}

\section{Setting a session-wide default}{
\code{bsparam()} will return a session-wide value for the BiocSingularParam object, analogous to \code{\link{bpparam}()}.
This defaults to a FastAutoParam object but can be modified by setting the \code{"BiocSingularParam.default"} global option to any BiocSingularParam object.
Any code that uses \code{bsparam()} will automatically use this new default without needing to explicitly pass the BiocSingularParam object to those functions.
}

\value{
Each constructor returns a BiocSingularParam subclass of the same type, containing the specified parameters.
}

\seealso{
\code{\link{runSVD}} for generic dispatch.

\code{\link{runExactSVD}}, \code{\link{runIrlbaSVD}} and \code{\link{runRandomSVD}} for specific methods.
}

\author{
Aaron Lun
}

\examples{
ExactParam()

IrlbaParam(tol=1e-8)

RandomParam(q=20)

# Modifying the default.
bsparam()
options(BiocSingularParam.default=IrlbaParam())
bsparam()
}
