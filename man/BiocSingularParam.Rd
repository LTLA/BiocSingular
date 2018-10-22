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

\title{BiocSingularParam classes}
\description{Classes for specifying the type of singular value decomposition (SVD) algorithm and associated parameters.}

\usage{
ExactParam(fold=5)

IrlbaParam(fold=5, extra.work=7, ...)

RandomParam(fold=5, ...)
    
\S4method{show}{BiocSingularParam}(object)

\S4method{show}{IrlbaParam}(object)

\S4method{show}{RandomParam}(object)
}

\arguments{
    \item{fold}{Numeric scalar, minimum fold-difference for cross-product calculation.}
    \item{extra.work}{Integer scalar, additional dimensionality of the workspace in \code{\link{runIrlba}}.}
    \item{...}{Additional arguments to pass to \code{\link{runIrlba}} or \code{\link{runRandomSVD}}.}
    \item{object}{An ExactParam object.}
}

\details{
The BiocSingularParam class is a virtual base class on which other parameter objects are built.
There are 3 concrete subclasses:
\describe{
    \item{}{\code{\link{ExactParam}}: exact SVD with \code{\link{runExactSVD}}.}
    \item{}{\code{\link{IrlbaParam}}: approximate SVD with \pkg{irlba} via \code{\link{runIrlba}}.}
    \item{}{\code{\link{RandomParam}}: approximate SVD with \pkg{rsvd} via \code{\link{runRandomSVD}}.}
}

These objects hold parameters specifying how each algorithm should be run on an arbitrary data set.
See the associated documentation pages for more details.
}

\section{Methods}{
\describe{
\item{\code{show(object)}:}{Display the class of a BiocSingularParam \code{object}.
Also shows the minimum fold-difference between dimensions required to perform the SVD on the cross-product rather than directly on the matrix.}
}
}

\seealso{
\code{\link{ExactParam}}, \code{\link{IrlbaParam}} and \code{\link{RandomParam}} for constructors.

\code{\link{runSVD}} for dispatch.
}

\author{
Aaron Lun
}

\examples{
ExactParam()

IrlbaParam(tol=1e-8)

RandomParam(q=20)
}
