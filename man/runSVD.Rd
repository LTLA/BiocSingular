\name{runSVD}
\alias{runSVD}
\alias{runSVD,missing-method}
\alias{runSVD,ExactParam-method}
\alias{runSVD,IrlbaParam-method}
\alias{runSVD,RandomParam-method}
\alias{runSVD,FastAutoParam-method}
\docType{methods}

\title{Run SVD}
\description{Perform a singular value decomposition on an input matrix with a specified algorithm.}

\usage{
runSVD(x, k, nu=k, nv=k, center=FALSE, scale=FALSE, 
    BPPARAM=SerialParam(), ..., BSPARAM=ExactParam())
}

\arguments{
\item{x}{A numeric matrix-like object to use in the SVD.}
\item{k}{Integer scalar specifying the number of singular values to return.}
\item{nu}{Integer scalar specifying the number of left singular vectors to return.}
\item{nv}{Integer scalar specifying the number of right singular vectors to return.}
\item{center}{Numeric vector, logical scalar or \code{NULL}, specifying values to subtract from each column of \code{x} - see \code{?"\link{BiocSingular-options}"}.}
\item{scale}{Numeric vector, logical scalar or \code{NULL}, specifying values to divide each column of \code{x} - see \code{?"\link{BiocSingular-options}"}.}
\item{BPPARAM}{A \link[BiocParallel]{BiocParallelParam} object specifying how parallelization should be performed.}
\item{...}{Further arguments to pass to specific methods.}
\item{BSPARAM}{A \linkS4class{BiocSingularParam} object specifying the type of algorithm to run.}
}

\details{
The class of \code{BSPARAM} will determine the algorithm that is used, see \code{?\linkS4class{BiocSingularParam}} for more details.
The default is to use an exact SVD via \code{\link{runExactSVD}}.
}

\value{
A list containing:
\itemize{
\item \code{d}, a numeric vector of the first \code{k} singular values.
\item \code{u}, a numeric matrix with \code{nrow(x)} rows and \code{nu} columns.
Each column contains a left singular vector.
\item \code{v}, a numeric matrix with \code{ncol(x)} rows and \code{nv} columns.
Each column contains a right singular vector.
}
}

\author{
Aaron Lun
}

\seealso{
\code{\link{runExactSVD}}, \code{\link{runIrlbaSVD}} and \code{\link{runRandomSVD}} for the specific functions.
}

\examples{
a <- matrix(rnorm(100000), ncol=20)

out.exact0 <- runSVD(a, k=4)
str(out.exact0)

out.exact <- runSVD(a, k=4, BSPARAM=ExactParam())
str(out.exact)

out.irlba <- runSVD(a, k=4, BSPARAM=IrlbaParam())
str(out.exact)

out.random <- runSVD(a, k=4, BSPARAM=RandomParam())
str(out.random)
}

