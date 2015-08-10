\name{hyper.like}
\alias{hyper.like}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The likelihood of a given set of data and an specified hyperplane
}
\description{
This is the mid-level likelihood solving function. Most users will not use this directly, but it is called by hyper.fit and hyper.plot2d/hyper.plot3d. Users can interact with the function directly, but it only takes arguments using the normal vector (coord.orth) orthogonal offset of the origin to the hyperplane (beta.orth) and orthogonal scatter to the hyperplane (scat.orth).
}
\usage{
hyper.like(parm, X, covarray, weights = 1, errorscale = 1, k.vec = FALSE, output = 'sum')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
A vector specifying the current paramters to compute the likelihood for. This should be a concatentation of 'normvec' coordinates and 'scat.orth' intrinsic scatter. e.g. for 3d data this would look like c(normvec1, normvec2, normvec3, scat.orth).
}
  \item{X}{
A position matrix with the N (number of data points) rows by d (number of dimensions) columns.
}
  \item{covarray}{
A dxdxN array containing the full covariance (d=dimensions, N=number of dxd matrices in the array stack). The 'makecovarray2d' and 'makecovarray3d' are convenience functions that make populating 2x2xN and 3x3xN arrays easier for a novice user.
}
  \item{weights}{
Vector of multiplicative weights for each row of the X data matrix. i.e. if this is 2 then it is equivalent to having two itentical data points with weights equal to 1. Should be either of length 1 (in which case elements are repeated as required) or the same length as the number of rows in the data matrix X.
}
 \item{errorscale}{
Value to multiplicatively re-scale the errors by (i.e. the covariance array become scaled by errorscale^2). This might be useful when trying to decide if the provided errors are too large.
}
\item{k.vec}{
A vector defining the direction of an exponential sampling distribution in the data. The length is the magnitude of the exponent, and it points in the direction of *increasing* density (see example below). Must be the same length as coord.orth, or FALSE.
}
  \item{output}{
If 'sum' then the output is the sum of the log likelihood. If 'val' then the output is the individual log likelihoods of the data provided. If 'sig' then the output represents the sigma tension of the data point with the current model, and can be thought of as -0.5chi^2 for each data point. 
}
}
\details{
hyper.convert is a convenience function that manipulates different parameterisations into the type required for the parm argument of hyper.like. See example below.
}
\value{
If output='sum' then the output is the sum of the log likelihood. If output='val' then the output is the individual log likelihoods of the data provided. If output='sig' then the output represents the sigma tension of the data point with the current model, and can be thought of as -0.5chi^2 for each data point. 
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
#Setup the initial data:

set.seed(650)
sampN=200
initscat=3
randatax=runif(sampN, -100,100)
randatay=rnorm(sampN, sd=initscat)
sx=runif(sampN, 0,10); sy=runif(sampN, 0,10)

mockvararray=makecovarray2d(sx, sy, corxy=0)

errxy={}
for(i in 1:sampN){
  rancovmat=ranrotcovmat2d(mockvararray[,,i])
  errxy=rbind(errxy, mvrnorm(1,mu=c(0,0), Sigma=rancovmat))
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

#Create our orthogonal vector. This does not need to be normalised to 1:

coord.orth=hyper.convert(coord=ang, in.coord.type = "theta", out.coord.type = "normvec")

#Feed this into the hyper.like function:

print(hyper.like(parm=coord.orth$parm, X, covarray, errorscale=1, output = "sum"))

#Comapre to a worse option:

print(hyper.like(parm=c(0.5, -1, 0, 4), X, covarray, errorscale=1, output = "sum"))

#As we can see, the paramters used to generate the data produce a higher likelihood.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }
\keyword{ hyper }% __ONLY ONE__ keyword per line
