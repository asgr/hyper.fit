\name{hyper.fit package}
\alias{hyper.fit-package}
\docType{package}
\title{
Hyperplane Fitting Package
}
\description{
Functions related to hyperplane fitting as presented in Robotham & Obreschkow 2015. Includes one main high level code for hyperplane fitting (hyper.fit). It also includes high level plotting and summary functions. In simple terms this allows the user to produce robust 1D linear fits for 2D x vs y type data, and robust 2D plane fits to 3D x vs y vs z type data. This hyperplane fitting works generically for any N-1 hyperplane model being fit to an N dimension data-set. All fits include intrinsic scatter in the generative hyperplane model.
}
\details{
\tabular{ll}{
Package: \tab hyper.fit\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 2016-08-01\cr
License: \tab GPL-3\cr
}
Most users will interact with this package through the high level hyper.fit and hyper.plot functions. The utility functions to construct covariance matrices and arrays will also be useful (see makecovmat2d and makecovarray2d).
}
\author{
Aaron Robotham and Danail Obreschkow

Maintainer: Aaron Robotham <aaron.robotham@uwa.edu.au>
}
\references{
Robotham, A.S.G., & Obreschkow, D., PASA, in press
}
\seealso{
  \code{\link{hyper.basic}}, \code{\link{hyper.convert}}, \code{\link{hyper.data}}, \code{\link{hyper.fit}}, \code{\link{hyper.plot}}, \code{\link{hyper.sigcor}}, \code{\link{hyper.summary}}
}
\keyword{ package }
\keyword{ fit }
\keyword{ hyper }
\keyword{ linear }
\keyword{ plane }
\keyword{ regression }

