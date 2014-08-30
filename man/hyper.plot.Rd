\name{hyper.plot}
\alias{hyper.plot2d}
\alias{hyper.plot3d}
\title{
A 2d and 3d likelihood diagnostic plot for optimal line fitting
}
\description{
These functions produce helpful 2d and 3d diagnostic plots for post hyper.fit analysis and for manual experimentation with parameter options. Error llipses and ellipsodis are added to the plots, with colouring scaled by 'sigma-tension' of the data points (where red is high tension). It also overplots the current line (2d) or plane (3d).
}
\usage{
hyper.plot2d(X, covarray, fitobj, parm.coord, parm.beta, parm.scat, coord.type='alpha', 
proj.type='orth', errorscale=1, clip=0.05, trans=1, ...)
hyper.plot3d(X, covarray, fitobj, parm.coord, parm.beta, parm.scat, coord.type='alpha', 
proj.type='orth', errorscale=1, clip=0.05, trans=1, ...)
%X,covarray,fitobj,parm.coord,parm.beta,parm.scat,coord.type='alpha',proj.type='orth',errorscale=1,clip=0.05,trans=0.5,...
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A position matrix with the N (numver of data points) rows by Dim (number of dimensions) columns.
}
  \item{covarray}{
A Dim x Dim x N array containing the full covariance. The 'makecovarray2d' and 'makecovarray3d' are convenience functions that make populating 2x2xN and 3x3xN arrays easier for a novice user.
}
  \item{fitobj}{
For simplicity the user can provide the direct output of hyper.fit to this argument. If this is not provided then parm.coord / parm.beta / parm.scat must all be specified.
}
  \item{parm.coord}{
Coord parameter/s to use. These are either angles that produce the vectors that form the N-dimensional slope or they are the vector elements directly.
}
  \item{parm.beta}{
Beta paramter to use. This is either specified as the distance from the origin to the hyperplane or as the intersection of the hyperplane on the fixed vertical axis of interest.
}
  \item{parm.scat}{
Intrinsic scatter parameter to use. This is either specified as the scatter orthogonal to the hyperplane or as the scatter along the vertical axis of interest.
}
  \item{coord.type}{
This specifies whether the parm.coord parameter is defined in terms of the unit vector of the line (alpha) or for the values of the angles that form the unit vector (theta).
}
  \item{proj.type}{
This specifies whether the parm.beta and the parm.scat are defined as orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis).
}
  \item{errorscale}{
Value to multiplicatively rescale the errors by (i.e. the covarince array become scaled by errorscale^2). This might be useful when trying to decide if the provided errors are too large.
}
  \item{clip}{
Quantile to colour clip the least likely data point. All data points below this value will be coloured red, and from there scaing progresses linearly via green to the most likely data point which will be blue.
}
  \item{trans}{
Transparency of the ellipses (hyper.plot2d) or ellipsoids (hyper.plot3d).
}
  \item{\dots}{
Arguments to pass to magplot (hyper.plot2d) or plot3d (hyper.plot3d).
}
}
\value{
The plotting functions also return the log-likelhood of the data points given the inputs.
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
hyper.plot2d(X, covarray, fitobj=fitline, trans=0.2)

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
newxyz=rotdata3d(randatax,randatay,randataz,theta=ang[1],dim='y')
newxyz=rotdata3d(newxyz[,1],newxyz[,2],newxyz[,3],theta=ang[2],dim='x')
mockplane=data.frame(x=newxyz[,1],y=newxyz[,2],z=newxyz[,3])

xerrang={};yerrang={};zerrang={}
corxyang={};corxzang={};coryzang={}
for(i in 1:sampN){
  newcovmatrot=rotcovmat(makecovmat3d(sx=sx[i], sy=sy[i], sz=sz[i], corxy=corxy[i],
  corxz=corxz[i],coryz=coryz[i]),theta=ang[1],dim='y')
  newcovmatrot=rotcovmat(newcovmatrot,theta=ang[2],dim='x')
  xerrang=c(xerrang,sqrt(newcovmatrot[1,1]))
  yerrang=c(yerrang,sqrt(newcovmatrot[2,2]))
  zerrang=c(zerrang,sqrt(newcovmatrot[3,3]))
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
covarray=makecovarray3d(mockplane$sx,mockplane$sy,mockplane$sz, mockplane$corxy,
mockplane$corxz, mockplane$coryz)
fitplane=hyper.fit(X=X, covarray=covarray, coord.type='theta',proj='orth')
print(fitplane)$parm
hyper.plot3d(X=X, covarray=covarray, fitplane, coord.type='theta',proj='orth')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ plot }
\keyword{ hyper }
\keyword{ linear }
\keyword{ plane }
\keyword{ regression }
