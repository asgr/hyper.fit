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
hyper.plot2d(X, covarray, vars, fitobj, parm.coord, parm.beta, parm.scat,
coord.type = 'alpha', proj.type = 'orth', errorscale = 1, doellipse = TRUE, clip = 0.05,
trans = 1, ...)
hyper.plot3d(X, covarray, vars, fitobj, parm.coord, parm.beta, parm.scat,
coord.type = 'alpha', proj.type = 'orth', errorscale = 1, doellipse = TRUE, clip = 0.05,
trans = 1, ...)
%X,covarray,fitobj,parm.coord,parm.beta,parm.scat,coord.type='alpha',proj.type='orth',errorscale=1,clip=0.05,trans=0.5,...
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A position matrix with the N (number of data points) rows by d (number of dimensions) columns.
}
  \item{covarray}{
A dxdxN array containing the full covariance (d=dimensions, N=number of dxd matrices in the array stack). The 'makecovarray2d' and 'makecovarray3d' are convenience functions that make populating 2x2xN and 3x3xN arrays easier for a novice user.
}
  \item{vars}{
A variance matrix with the N (numver of data points) rows by Dim (number of dimensions) columns. In effect this is the diagonal elements of the 'covarray' array where all other terms are zero. If 'covarray' is also provided that is used instead.
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
  \item{doellipse}{
Should 2d/3d error ellipses be drawn on the plot? This should only be TRUE if data errors are actually provided, else it should be FALSE.  
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
hyper.plot2d(X=X, covarray=covarray, fitobj=fitline, trans=0.2, asp=1)

################### A 3D example with fitting a plane ###################

\dontrun{

#Setup the initial data:

set.seed(650)
sampN=200
initscat=3
randatax=runif(sampN, -100, 100)
randatay=runif(sampN, -100, 100)
randataz=rnorm(sampN, sd=initscat)
sx=runif(sampN, 0, 5); sy=runif(sampN,0,5); sz=runif(sampN, 0, 5)

mockvararray=makecovarray3d(sx, sy, sz, corxy=0, corxz=0, coryz=0)

errxyz={}
for(i in 1:sampN){
  rancovmat=ranrotcovmat3d(mockvararray[,,i])
  errxyz=rbind(errxyz, mvrnorm(1, mu=c(0, 0, 0), Sigma=rancovmat))
  mockvararray[,,i]=rancovmat
  }
randatax=randatax+errxyz[,1]
randatay=randatay+errxyz[,2]
randataz=randataz+errxyz[,3]
sx=sqrt(mockvararray[1,1,]); sy=sqrt(mockvararray[2,2,]); sz=sqrt(mockvararray[3,3,])
corxy=mockvararray[1,2,]/(sx*sy); corxz=mockvararray[1,3,]/(sx*sz)
coryz=mockvararray[2,3,]/(sy*sz)

#Rotate the data to an arbitrary angle theta/phi:
desiredxtozang=10
desiredytozang=40
ang=c(desiredxtozang*cos(desiredytozang*pi/180), desiredytozang)
newxyz=rotdata3d(randatax, randatay, randataz, theta=ang[1], dim='y')
newxyz=rotdata3d(newxyz[,1], newxyz[,2], newxyz[,3], theta=ang[2], dim='x')
mockplane=data.frame(x=newxyz[,1], y=newxyz[,2], z=newxyz[,3])

xerrang={}; yerrang={}; zerrang={}
corxyang={}; corxzang={}; coryzang={}
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

X=cbind(mockplane$x, mockplane$y, mockplane$z)
covarray=makecovarray3d(mockplane$sx, mockplane$sy, mockplane$sz, mockplane$corxy,
mockplane$corxz, mockplane$coryz)
fitplane=hyper.fit(X=X, covarray=covarray, coord.type='theta', proj='orth')
hyper.plot3d(X=X, covarray=covarray, fitobj=fitplane)
}

################### Example using the data from Hogg 2010 ###################

#Example using the data from Hogg 2010: http://arxiv.org/pdf/1008.4686v1.pdf

#Full data

hogg=read.table(system.file('data/hogg.tab', package='hyper.fit'), header=TRUE)
fithogg=hyper.fit(X=cbind(hogg$x, hogg$y), covarray=makecovarray2d(hogg$sx, hogg$sy,
hogg$corxy), coord.type='theta', proj.type='orth')
hyper.plot2d(X=cbind(hogg$x, hogg$y), covarray=makecovarray2d(hogg$sx, hogg$sy,
hogg$corxy), fitobj=fithogg, trans=0.2, xlim=c(0, 300), ylim=c(0, 700))

#We now do exercise 17 of Hogg 2010 using trimmed data, where we remove the high tension
#data point 3 (which we can see as the reddest point in the above plot:

hogg=read.table(system.file('data/hogg.tab', package='hyper.fit'), header=TRUE)
hoggtrim=hogg[-3,]
fithoggtrim=hyper.fit(X=cbind(hoggtrim$x, hoggtrim$y), covarray=makecovarray2d(hoggtrim$sx,
hoggtrim$sy, hoggtrim$corxy), coord.type='theta', proj.type='orth', algo.func='LA')
hyper.plot2d(X=cbind(hoggtrim$x, hoggtrim$y), covarray=makecovarray2d(hoggtrim$sx,
hoggtrim$sy, hoggtrim$corxy), fitobj=fithoggtrim, trans=0.2, xlim=c(0, 300), ylim=c(0, 700))

#We can compare this against the previous fit with:
hyper.plot2d(cbind(hoggtrim$x, hoggtrim$y), covarray=makecovarray2d(hoggtrim$sx,
hoggtrim$sy, hoggtrim$corxy), fitobj=fithogg, trans=0.2, xlim=c(0, 300), ylim=c(0, 700))

################### Example using 'real' data with intrinsic scatter ###################

intrin=read.table(system.file('data/intrin.tab', package='hyper.fit'), header=TRUE)
fitintrin=hyper.fit(X=cbind(intrin$x, intrin$y), covarray=makecovarray2d(intrin$sx,
intrin$sy, intrin$corxy), coord.type='theta', proj.type='orth')
hyper.plot2d(cbind(intrin$x, intrin$y), covarray=makecovarray2d(intrin$sx,
intrin$sy, intrin$corxy), fitobj=fitintrin, trans=0.1, pch='.', asp=1)

################### Example using flaring trumpet data ###################

trumpet=read.table(system.file('data/trumpet.tab', package='hyper.fit'), header=TRUE)
fittrumpet=hyper.fit(X=cbind(trumpet$x, trumpet$y), covarray=makecovarray2d(trumpet$sx,
trumpet$sy,trumpet$corxy), coord.type='theta')
hyper.plot2d(cbind(trumpet$x, trumpet$y), covarray=makecovarray2d(trumpet$sx,
trumpet$sy, trumpet$corxy), fitobj=fittrumpet, trans=0.1, pch='.', asp=1)
#If you look at the ?hyper.fit example we find that zero intrinsic scatter is actually
#preferred, but we don't see this in the above plot.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ plot }
\keyword{ hyper }
\keyword{ linear }
\keyword{ plane }
\keyword{ regression }
