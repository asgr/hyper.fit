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
hyper.like(X, covarray, coord.orth, beta.orth = 0, scat.orth = 1, errorscale=1,
output = "sum")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A position matrix with the N (numver of data points) rows by Dim (number of dimensions) columns.
}
  \item{covarray}{
A Dim x Dim x N array containing the full covariance. The 'makecovarray2d' and 'makecovarray3d' are convenience functions that make populating 2x2xN and 3x3xN arrays easier for a novice user.
}
  \item{coord.orth}{
The vector that defines the normal to the hyper-plane. This does not need to be normalised.
}
  \item{beta.orth}{
The orthogonal distance between the origina and the hyperplane.
}
  \item{scat.orth}{
The intrinsic scatter orthogonal to the hyperplane.
}
 \item{errorscale}{
Value to multiplicatively rescale the errors by (i.e. the covarince array become scaled by errorscale^2). This might be useful when trying to decide if the provided errors are too large.
}
  \item{output}{
If 'sum' then the output is the sum of the log likelihood. If 'val' then the output is the individual log likelihoods of the data provided. If 'sig' then the output represents the sigma tension of the data point with the current model, and can be thought of as -0.5chi^2 for each data point. 
}
}
\details{
'coord.orth' is slightly different to the raw ouput of the fit or parm.coord input in that the reference axis must explitly have a value of -1 in the relevant position, i.e. y=2x becomes coord.orth=c(2,-1) and z=2x+3y becomes coord.orth=c(2,3,-1). If the reference axis is not the last then care must be taken, e.g. y=2x+3z becomes coord.orth=c(2,-1,3) since the reference axis y is the second element in the normal vector. Clearly dividing through by -3 and rearranging takes us back to the more familiar z=y/3-2x/3 and coord.orth=c(-2/3,1/3,-1), which we get to just by dividing through by -3.
}
\value{
If output='sum' then the output is the sum of the log likelihood. If output='val' then the output is the individual log likelihoods of the data provided. If output='sig' then the output represents the sigma tension of the data point with the current model, and can be thought of as -0.5chi^2 for each data point. 
}
\references{
Robotham & Obreschkow 2014
}
\author{
Aaron Robotham and Danail Obreschkow
}
\seealso{
\code{\link{hyper.fit}}, \code{\link{hyper.plot2d}}, \code{\link{hyper.plot3d}}, \code{\link{makecovmat2d}}, \code{\link{makecovmat3d}}
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

coord.orth=c(coord.convert(ang, in.coord.type = "theta", out.coord.type = "alpha"),-1)

#Feed this into the hyepr.like function:

print(hyper.like(X, covarray, coord.orth=coord.orth, beta.orth = 0, scat.orth = 3,
errorscale=1, output = "sum"))

#Comapre to a worse option:

print(hyper.like(X, covarray, coord.orth=c(0.5,-1), beta.orth = 0, scat.orth = 4,
errorscale=1, output = "sum"))

#As we can see, the paramters used to generate the data produce a higher likelihood.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }
\keyword{ hyper }% __ONLY ONE__ keyword per line
