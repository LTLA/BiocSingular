\name{runSVD}
\alias{runSVD}
\alias{runSVD,missing-method}
\alias{runSVD,NULL-method}
\alias{runSVD,ExactParam-method}
\alias{runSVD,IrlbaParam-method}
\alias{runSVD,RandomParam-method}
\docType{methods}

\title{Run SVD}
\description{Perform a singular value decomposition on an input matrix with a specified algorithm.}

\usage{
runSVD(..., BSPARAM=NULL)
}

\arguments{
\item{...}{Further arguments to pass to individual functions.
This can contain:
    \itemize{
        \item{\code{x}: A numeric matrix-like object to use in the SVD.}
        \item{\code{k}: Integer scalar specifying the number of singular values to return.}
        \item{\code{nu}: Integer scalar specifying the number of left singular vectors to return.}
        \item{\code{nv}: Integer scalar specifying the number of right singular vectors to return.}
        \item{\code{center}: Numeric vector of length equal to \code{ncol(x)}, specifying values to subtract from each column of \code{x}.}
        \item{\code{scale}: Numeric vector of length equal to \code{ncol(x)}, specifying values to divide each column of \code{x}.}
        \item{\code{BPPARAM}: A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.}
    }
Of these, only \code{x} is strictly mandatory, though specifying \code{k}, \code{nu} or \code{nv} is often necessary in practical use. 
}
\item{BSPARAM}{A \linkS4class{BiocSingularParam} object specifying the type of algorithm to run.}
}

\details{
The class of \code{BSPARAM} will determine the algorithm that is used, see \code{?\linkS4class{BiocSingularParam}} for more details.

If \code{BSPARAM=NULL} or is missing, the method will default to an exact SVD via \code{\link{runExactSVD}}.
}

\value{
A list containing:
\itemize{
\item \code{d}, a numeric vector of the first \code{k} singular values.
\item \code{u}, a numeric matrix with \code{nrow(x)} rows and \code{nu} columns.
Each column contains a left singular vector.
\item \code{u}, a numeric matrix with \code{ncol(x)} rows and \code{nv} columns.
Each column contains a right singular vector.
}
}

\author{
Aaron Lun
}

\seealso{
\code{\link{runExactSVD}}, \code{\link{runIrlba}} and \code{\link{runRandomSVD}} for the specific functions.
}

\examples{
a <- matrix(rnorm(100000), ncol=20)

out.exact0 <- runSVD(a)
str(out.exact0)

out.exact <- runSVD(a, BSPARAM=ExactParam())
str(out.exact)

out.irlba <- runSVD(a, BSPARAM=IrlbaParam())
str(out.exact)

out.random <- runSVD(a, BSPARAM=RandomParam())
str(out.random)
}

