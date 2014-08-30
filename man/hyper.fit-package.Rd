\name{hyper.fit-package}
\alias{hyper.fit-package}
\docType{package}
\title{
Hyperplane fitting with intrinsic scatter
}
\description{
This package includes two main high level codes for hyperplane fitting (hyper.fit) and visualising (hyper.plot2d / hyper.plot3d). In simple terms this allows the user to produce robust 1D linear fits for 2D x vs y type data, and robust 2D plane fits to 3D x vs y vs z type data. This hyperplane fitting works generically for any N-1 hyperplane model being fit to a N dimension dataset. All fits include intrinsic scatter in the generative model orthogonal to the hyperplane.
}
\details{
\tabular{ll}{
Package: \tab hyper.fit\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-08-29\cr
License: \tab GPL-2\cr
}
Most users will interact with this package through the high level hyper.fit and hyper.plot functions. The utility functions to construct covariance matrices and arrays will also be useful, see makecovmat2d and makecovarray2d.
}
\author{
Aaron Robotham and Danail Obreschkow

Maintainer: Aaron Robotham <aaron.robotham@uwa.edu.au>
}
\references{
Robotham & Obreschkow 2014
}
\keyword{ package }
\keyword{ fit }
\keyword{ hyper }
\keyword{ linear }
\keyword{ plane }
\keyword{ regression }
\seealso{
\code{\link{hyper.fit}}, \code{\link{hyper.plot2d}}, \code{\link{hyper.plot3d}}, \code{\link{makecovmat2d}}, \code{\link{makecovmat3d}}
}

