\name{Sigma Correct}
\alias{hyper.sigcor}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to convert from biased sample sigma to unbiased population sigma
}
\description{
  Calculates the required corrections for transforming the biased estimate of sigma to an unbiased estimate, and for transforming the sample expectation to the population expectation.
  
}
\usage{
  hyper.sigcor(N,parmDF)
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
  A vector of length 3. The first element is the correction from the biased to unbiased estimate for sigma. The second element is the correction from the sample to population estimate for sigma. The third is the combination of the previous two (i.e. the total correction the user will typically want to apply to their data).
  
}
\references{
  Robotham, A.S.G., & Obreschkow, D., 2014
}
\author{
  Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.basic}}, \code{\link{hyper.convert}}, \code{\link{hyper.data}}, \code{\link{hyper.fit}}, \code{\link{hyper.plot}}, \code{\link{hyper.sigcor}}, \code{\link{hyper.summary}}
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
