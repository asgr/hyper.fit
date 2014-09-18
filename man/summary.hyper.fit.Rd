\name{summary.hyper.fit}
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
  Robotham, A.S.G., & Obreschkow, D., 2014
}
\author{
  Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.fit}}, \code{\link{hyper.plot2d}}, \code{\link{hyper.plot3d}}, \code{\link{makecovmat2d}}, \code{\link{makecovmat3d}}
}

\examples{
################### A 2D example with fitting a line ###################

#Setup the initial data:

set.seed(650)
sampN=200
initscat=3
randatax=runif(sampN, -100, 100)
randatay=rnorm(sampN, sd=initscat)
sx=runif(sampN, 0, 10); sy=runif(sampN, 0, 10)

mockvararray=makecovarray2d(sx, sy, corxy=0)

errxy={}
for(i in 1:sampN){
  rancovmat=ranrotcovmat2d(mockvararray[,,i])
  errxy=rbind(errxy, mvrnorm(1, mu=c(0, 0), Sigma=rancovmat))
  mockvararray[,,i]=rancovmat
  }
randatax=randatax+errxy[,1]
randatay=randatay+errxy[,2]

#Rotate the data to an arbitrary angle theta:

ang=30
mock=rotdata2d(randatax, randatay, theta=ang)
xerrang={}; yerrang={}; corxyang={}
for(i in 1:sampN){
  covmatrot=rotcovmat(mockvararray[,,i], theta=ang)
  xerrang=c(xerrang, sqrt(covmatrot[1,1])); yerrang=c(yerrang, sqrt(covmatrot[2,2]))
  corxyang=c(corxyang, covmatrot[1,2]/(xerrang[i]*yerrang[i]))
}
corxyang[xerrang==0 & yerrang==0]=0
mock=data.frame(x=mock[,1], y=mock[,2], sx=xerrang, sy=yerrang, corxy=corxyang)

#Do the fit:

X=cbind(mock$x, mock$y)
covarray=makecovarray2d(mock$sx, mock$sy, mock$corxy)
fitline=hyper.fit(X=X, covarray=covarray, coord.type='theta')
summary(fitline)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ summary }% __ONLY ONE__ keyword per line
