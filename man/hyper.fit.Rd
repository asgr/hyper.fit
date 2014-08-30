\name{hyper.fit}
\alias{hyper.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Top level function that attempts to fit a hyperplane to provided data.
}
\description{
Top level line fitting function that uses downhill searches (optim/LaplaceApproximation) or MCMC (LaplacesDemon) to search out the best fitting parameters for a hyperplane (minimum a 1D line for 2D data), including the intrinsic scatter as part of the fit.
}
\usage{
hyper.fit(X, covarray, vars, parm.coord, parm.beta, parm.scat, vert.axis, itermax = 1e4,
coord.type = 'alpha', proj.type = 'orth', algo.func = 'optim', algo.method = 'default',
Specs = list(alpha.star=0.44), doerrorscale=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A position matrix with the N (numver of data points) rows by Dim (number of dimensions) columns.
}
  \item{covarray}{
A Dim x Dim x N array containing the full covariance. The 'makecovarray2d' and 'makecovarray3d' are convenience functions that make populating 2x2xN and 3x3xN arrays easier for a novice user.
}
  \item{vars}{
A variance matrix with the N (numver of data points) rows by Dim (number of dimensions) columns. In effect this is the diagonal elements of the 'covarray' array where all other terms are zero. If 'covarray' is also provided that is used instead.
}
  \item{parm.coord}{
Vector of initial coord paramters. These are either angles that produce the vectors that form the N-dimensional slope or they are the vector elements directly.
}
  \item{parm.beta}{
Initial value of beta. This is either specified as the distance from the origin to the hyperplane or as the intersection of the hyperplane on the fixed vertical axis of interest.
}
  \item{parm.scat}{
Initial value of the intrinsic scatter. This is either specified as the scatter orthogonal to the hyperplane or as the scatter along the vertical axis of interest.
}
  \item{vert.axis}{
Which axis should the plane equation be formulated for. This must be a number which specifies the column of position matrix 'X' to be defined by the plane.
}
  \item{itermax}{
The maximum iterations to use for either the LaplaceApproximation function or LaplacesDemon function.
}
  \item{coord.type}{
This specifies whether the fit should be done in terms of the unit vector of the line (alpha) or by the values of the angles that form the unit vector (theta).
}
  \item{proj.type}{
This specifies whether the fit offset beta and the intrinsic scatter should be defined orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis).
}
  \item{algo.func}{
If 'optim' will optimise using the R base optim function. If 'LA' will optimise using the LaplaceApproximation function. If 'LD' will optimise using the LaplacesDemon function. For both 'LA' and 'LA' the 'LaplacesDemon' package will be used (see http://www.bayesian-inference.com/software).
}
  \item{algo.method}{
Specifies 'method' when using 'optim'. Specifies 'Method' when using 'LaplaceApproximation'. Specifies 'Algorithm' when using 'LaplacesDemon'.
}
  \item{Specs}{
Specs to pass to the 'LaplacesDemon' function. Default options are for the default CHARM algorithm.
}
  \item{doerrorscale}{
If FALSE then the provided covariance are treated as it. If TRUE then the likelihood function is also allowed to rescale all errors my a uniform amount.
}
}
\details{
Setting errorscale to TRUE allows for hyper stable solutions when errors are erroneously inflated, e.g. when even with intrinsic scatter equal to one the data is more clustered around the optimal likelihood plane then you should expect given the provided covariance matrix. See Examples below for a 2D scenario where this is helpful.

algo.func='LD' also provides the chance that the true generative model has an intrinsic scatter of zero. The other available functions (optim and LA) can find exact solution equal to zero since they are strictly mode finding (i.e. maximum likelihood) routines. LD is MCMC based, so naturally returns a mean/expectation which *must* have a finite positive value for the intrinsic scatter (it's not allowed to travel below zero).
}
\value{

The function returns a multi-component list containing:

\item{parm}{
  The main paramter fit outputs specified as set by the coord.type and proj.type options.
}
\item{parm.vert.axis}{
  The main paramter fit outputs specified strictly along the defined vert.axis (for both the intrinsic scatter and the offset).
}
\item{covar}{
  The covariance matrix for parm. Only for algo.func='optim' and algo.func='LA'.
}
\item{covar.vert.axis}{
  The covariance matrix for parm.vert.axis, specified strictly along the defined vert.axis (for both the intrinsic scatter and the offset). Only for algo.func='optim' and algo.func='LA'.
}
\item{fit}{
  The direct output of the specified algo.func. So either the natural return from optim, LaplaceApproximation or LaplacesDemon.
}
\item{zeroscatprob}{
  The fraction of samples for the intrinsic scatter which are at *exactly* zero, which provides a guideline probability for the intrinsic scatter being truly zero, rather than the expectation which will always be a finite amount above zero. Only for algo.func='LD'.
}
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
randatax=runif(sampN,-100,100)
randatay=rnorm(sampN,sd=initscat)
sx=runif(sampN,0,10);sy=runif(sampN,0,10)

mockvararray=makecovarray2d(sx,sy,corxy=0)

errxy={}
for(i in 1:sampN){
  rancovmat=ranrotcovmat2d(mockvararray[,,i])
  errxy=rbind(errxy,mvrnorm(1,mu=c(0,0),Sigma=rancovmat))
  mockvararray[,,i]=rancovmat
  }
randatax=randatax+errxy[,1]
randatay=randatay+errxy[,2]

#Rotate the data to an arbitrary angle theta:

ang=30
mock=rotdata2d(randatax,randatay,theta=ang)
xerrang={};yerrang={};corxyang={}
for(i in 1:sampN){
  covmatrot=rotcovmat(mockvararray[,,i],theta=ang)
  xerrang=c(xerrang,sqrt(covmatrot[1,1]));yerrang=c(yerrang,sqrt(covmatrot[2,2]))
  corxyang=c(corxyang,covmatrot[1,2]/(xerrang[i]*yerrang[i]))
}
corxyang[xerrang==0 & yerrang==0]=0
mock=data.frame(x=mock[,1],y=mock[,2],sx=xerrang,sy=yerrang,corxy=corxyang)

#Do the fit:

X=cbind(mock$x,mock$y)
covarray=makecovarray2d(mock$sx,mock$sy,mock$corxy)
fitline=hyper.fit(X, covarray=covarray, coord.type='theta')
print(fitline$parm)

#Now a 3D example with fitting a plane:

set.seed(650)
sampN=200
initscat=3
randatax=runif(sampN,-100,100)
randatay=runif(sampN,-100,100)
randataz=rnorm(sampN,sd=initscat)
sx=runif(sampN,0,5);sy=runif(sampN,0,5);sz=runif(sampN,0,5)

mockvararray=makecovarray3d(sx,sy,sz,corxy=0,corxz=0,coryz=0)

errxyz={}
for(i in 1:sampN){
  rancovmat=ranrotcovmat3d(mockvararray[,,i])
  errxyz=rbind(errxyz,mvrnorm(1,mu=c(0,0,0),Sigma=rancovmat))
  mockvararray[,,i]=rancovmat
  }
randatax=randatax+errxyz[,1]
randatay=randatay+errxyz[,2]
randataz=randataz+errxyz[,3]
sx=sqrt(mockvararray[1,1,]);sy=sqrt(mockvararray[2,2,]);sz=sqrt(mockvararray[3,3,])
corxy=mockvararray[1,2,]/(sx*sy); corxz=mockvararray[1,3,]/(sx*sz)
coryz=mockvararray[2,3,]/(sy*sz)

#Rotate the data to an arbitrary angle theta/phi:
desiredxtozang=10
desiredytozang=40
ang=c(desiredxtozang*cos(desiredytozang*pi/180),desiredytozang)
newxyz=rotdata3d(randatax, randatay, randataz, theta=ang[1], dim='y')
newxyz=rotdata3d(newxyz[,1], newxyz[,2], newxyz[,3], theta=ang[2], dim='x')
mockplane=data.frame(x=newxyz[,1], y=newxyz[,2], z=newxyz[,3])

xerrang={};yerrang={};zerrang={}
corxyang={};corxzang={};coryzang={}
for(i in 1:sampN){
  newcovmatrot=rotcovmat(makecovmat3d(sx=sx[i], sy=sy[i], sz=sz[i], corxy=corxy[i],
  corxz=corxz[i], coryz=coryz[i]), theta=ang[1], dim='y')
  newcovmatrot=rotcovmat(newcovmatrot, theta=ang[2], dim='x')
  xerrang=c(xerrang, sqrt(newcovmatrot[1,1]))
  yerrang=c(yerrang, sqrt(newcovmatrot[2,2]))
  zerrang=c(zerrang, sqrt(newcovmatrot[3,3]))
  corxyang=c(corxyang, newcovmatrot[1,2]/(xerrang[i]*yerrang[i]))
  corxzang=c(corxzang, newcovmatrot[1,3]/(xerrang[i]*zerrang[i]))
  coryzang=c(coryzang, newcovmatrot[2,3]/(yerrang[i]*zerrang[i]))
}
corxyang[xerrang==0 & yerrang==0]=0
corxzang[xerrang==0 & zerrang==0]=0
coryzang[yerrang==0 & zerrang==0]=0
mockplane=data.frame(x=mockplane$x, y=mockplane$y, z=mockplane$z, sx=xerrang, sy=yerrang,
sz=zerrang, corxy=corxyang, corxz=corxzang, coryz=coryzang)

X=cbind(mockplane$x,mockplane$y,mockplane$z)
covarray=makecovarray3d(mockplane$sx, mockplane$sy, mockplane$sz, mockplane$corxy,
mockplane$corxz, mockplane$coryz)
fitplane=hyper.fit(X=X, covarray=covarray, coord.type='theta',proj='orth')
print(fitplane)$parm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ hyper }
\keyword{ linear }
\keyword{ plane }
\keyword{ regression }
