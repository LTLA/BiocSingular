\name{ResidualMatrix}

\alias{ResidualMatrixSeed}
\alias{ResidualMatrix}

\docType{class}

\title{The ResidualMatrix class}
\description{
This class is deprecated, see the exact same class in the \pkg{ResidualMatrix} package.
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

\author{
Aaron Lun
}
