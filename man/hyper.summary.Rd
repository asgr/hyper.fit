\name{hyper.summary}
\alias{hyper.summary}
\alias{summary.hyper.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary function for hyper.fit object
}
\description{
Prints out basic summary information for hyper.fit objects output by the hyper.fit function.
}
\usage{
\method{summary}{hyper.fit}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class hyper.fit. This is the only structure that needs to be provided when executing summary(fitobj) class dependent plotting, which will use the summary.hyper.fit function.  
}
  \item{\dots}{
  Arguments passed to summary function.
}
}
\details{
  Outputs basic summary of the hyper.fit output.
}
\value{
  Prints various summary outputs.
}
\references{
  Robotham, A.S.G., & Obreschkow, D., PASA, in press
}
\author{
  Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.basic}}, \code{\link{hyper.convert}}, \code{\link{hyper.data}}, \code{\link{hyper.fit}}, \code{\link{hyper.plot}}, \code{\link{hyper.sigcor}}, \code{\link{hyper.summary}}
}

\examples{
#### Example using 6dFGS Fundamental Plane data ####

FP6dFGS=read.table(system.file('data/FP6dFGS.tab', package='hyper.fit'), header=TRUE)
fitFP6dFGSw=hyper.fit(FP6dFGS[,c('logIe_J', 'logsigma', 'logRe_J')],
vars=FP6dFGS[,c('logIe_J_err', 'logsigma_err', 'logRe_J_err')]^2, weights=FP6dFGS[,'weights'],
coord.type='alpha', scat.type='vert.axis')
summary(fitFP6dFGSw)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ summary }% __ONLY ONE__ keyword per line
