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

\docType{class}
\alias{show,BiocSingularParam-method}

\alias{bsfold}
\alias{bsdeferred}

\title{BiocSingularParam classes}
\description{Classes for specifying the type of singular value decomposition (SVD) algorithm and associated parameters.}

\usage{
ExactParam(deferred=FALSE, fold=5)

IrlbaParam(deferred=FALSE, fold=5, extra.work=7, ...)

RandomParam(deferred=FALSE, fold=5, ...)
    
bsfold(object)

bsdeferred(object)
}

\arguments{
    \item{deferred}{Logical scalar indicating whether centering/scaling should be deferred, see \code{?"\link{BiocSingular-options}"}.}
    \item{fold}{Numeric scalar specifying the minimum fold-difference for cross-product calculation, see \code{?"\link{BiocSingular-options}"}.}
    \item{extra.work}{Integer scalar, additional dimensionality of the workspace in \code{\link{runIrlbaSVD}}.}
    \item{...}{Additional arguments to pass to \code{\link{runIrlbaSVD}} or \code{\link{runRandomSVD}}.}
    \item{object}{A BiocSingularParam object.}
}

\details{
The BiocSingularParam class is a virtual base class on which other parameter objects are built.
There are 3 concrete subclasses:
\describe{
    \item{}{\code{\link{ExactParam}}: exact SVD with \code{\link{runExactSVD}}.}
    \item{}{\code{\link{IrlbaParam}}: approximate SVD with \pkg{irlba} via \code{\link{runIrlbaSVD}}.}
    \item{}{\code{\link{RandomParam}}: approximate SVD with \pkg{rsvd} via \code{\link{runRandomSVD}}.}
}

These objects hold parameters specifying how each algorithm should be run on an arbitrary data set.
See the associated documentation pages for more details.
}

\section{Methods}{
\describe{
\item{\code{show(object)}:}{Display the class of a BiocSingularParam \code{object}, and a summary of the set parameters.}
\item{\code{bsfold(object)}:}{Return a numeric scalar specifying the fold-difference for cross-product calculation.}
\item{\code{bsdeferred(object)}:}{Return a logical scalar indicating whether centering and scaling should be deferred.}
}
}

\value{
The constructors return a BiocSingularParam subclass of the same type, containing the specified parameters.
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
}
