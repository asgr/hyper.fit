\name{Sigma Corrections}
\alias{sigcor}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to convert from biased sample sigma to unbiased population sigma
}
\description{
  sigcor: Calculates the required corrections for transforming the biased estimate of sigma to an unbiased estimate, and for transforming the sample expectation to the population expectation.
  
}
\usage{
  sigcor(N,parmDF)
}
\arguments{
  \item{N}{
    The number of data points in the sample where sigma is being estimated.  
  }
  \item{parmDF}{
    The number of degrees of freedom for the parameters. In the case of fitting a 1d Gaussian to data via maximum likelihood (equivalent the measuring the standard deviation of the data) parmDF=1, when fitting a 1d line with intrinsic scatter parmDF=2 and when fitting a plane to 3d data with intrinsic scatter parmDF=3.
  }
}
\value{
  A vector of length 3. The first element is the correction from the biased to unbiased estimate for sigma. The second element is the correction from the sample to population estimate for sigmae. The third is the combination of the previous two (i.e. the total correction the suer will typically want to apply to their data).
  
}
\references{
  Robotham, A.S.G., & Obreschkow, D., 2014
}
\author{
  Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.fit}}, \code{\link{hyper.plot2d}}, \code{\link{hyper.plot3d}}, \code{\link{makecovmat2d}}, \code{\link{makecovmat3d}}
}
\examples{
  #None yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ utility }
\keyword{ sigma }
\keyword{ bias }
\keyword{ sample }
\keyword{ population }
