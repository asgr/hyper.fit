\name{hyper.plot}
\alias{hyper.plot}
\alias{hyper.plot2d}
\alias{hyper.plot3d}
\alias{plot.hyper.fit}
\title{
A 2d and 3d likelihood diagnostic plot for optimal line fitting
}
\description{
These functions produce helpful 2d and 3d diagnostic plots for post hyper.fit analysis and for manual experimentation with parameter options. Error llipses and ellipsodis are added to the plots, with colouring scaled by 'sigma-tension' of the data points (where red is high tension). It also overplots the current line (2d) or plane (3d). If the data is either 2d/3d then a simple interface to the relevant plot.hyper function is provided by the hyper.fit class dependent function plot.hyper.fit, where the user only has to execute plot(fitoutput).
}
\usage{
\method{plot}{hyper.fit}(x, ...)

hyper.plot2d(X, covarray, vars, fitobj, parm.coord, parm.beta, parm.scat,
parm.errorscale = 1, vert.axis, weights, k.vec, coord.type = 'alpha', scat.type = 'orth',
doellipse = TRUE, sigscale=c(0,4), trans=1, dobar=FALSE, position='topright', ...)

hyper.plot3d(X, covarray, vars, fitobj, parm.coord, parm.beta, parm.scat,
parm.errorscale = 1, vert.axis, weights, k.vec, coord.type = 'alpha', scat.type = 'orth',
doellipse = TRUE, sigscale=c(0,4),trans=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Argument for the class dependent plot.hyper.fit function. An object of class hyper.fit. This is the only structure that needs to be provided when executing plot(fitobj) class dependent plotting, which will use the plot.hyper.fit function.  
}
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
For simplicity the user can provide the direct output of hyper.fit to this argument, which sets the hyperplane parameter values to those found during the fitting process. If this is not provided then parm.coord / parm.beta / parm.scat must all be specified.
}
  \item{parm.coord}{
Vector of initial coord paramters. These are either angles that produce the vectors that predict the vert.axis dimension (coord.type='theta'), the gradients of these (coord.type='alpha') or they are the unit vector normal to the hyperplane (coord.type='unitvec').
}
  \item{parm.beta}{
Initial value of beta. This is either specified as the absolute distance from the origin to the hyperplane or as the intersection of the hyperplane on the vert.axis dimension being predicted.
}
  \item{parm.scat}{
Initial value of the intrinsic scatter. This is either specified as the scatter orthogonal to the hyperplane or as the scatter along the vert.axis dimension being predicted.
}
  \item{parm.errorscale}{
Value to multiplicatively rescale the errors by (i.e. the covarince array becomes scaled by errorscale^2). This might be useful when trying to decide if the provided errors are too large. Default is 1, and therefore it does not need to be specified explicitly.
}
  \item{vert.axis}{
Which axis should the plane equation be formulated for. This must be a number which specifies the column of position matrix 'X' to be defined by the hyperplane. If missing, then the projection dimension is assumed to be the last column of the 'X' matrix.
}
  \item{weights}{
Vector of multiplicative weights for each row of the X data matrix. i.e. if this is 2 then it is equivalent to having two itentical data points with weights equal to 1. Should be either of length 1 (in which case elements are repeated as required) or the same length as the number of rows in the data matrix X.
}
  \item{k.vec}{
A vector defining the direction of an exponential sampling distribution in the data. The length is the scaling 'a' of the generative exponent (i.e., exp(a*x)), and it points in the direction of *increasing* density (see example below). If provided, k.vec must be the same length as the dimensions of the data. k.vec has the most noticeable effect on the beta offset parameter. It is correcting for the phenomenom Astronomers often call Eddington bias.
}
  \item{coord.type}{
This specifies whether the parm.coord parameter is defined in terms of the unit vector of the line (alpha) or for the values of the angles that form the unit vector (theta).
}
  \item{scat.type}{
This specifies whether the parm.beta and the parm.scat are defined as orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis).
}
  \item{doellipse}{
Should 2d/3d error ellipses be drawn on the plot? This should only be TRUE if data errors are actually provided, else it should be FALSE.  
}
  \item{sigscale}{
Vector of length 2 specifying the lower and upper limits for the linear blue->green->red mapping used to colour the data. The default will map to 0->2->4 sigma offset tension. 
}
  \item{trans}{
Transparency of the ellipses (hyper.plot2d) or ellipsoids (hyper.plot3d).
}
  \item{dobar}{
Logical specifying whether a magbar colour scale is added to the 2D plot. Only available for hyper.plot2d.
}
  \item{position}{
If dobar=TRUE, then position specifies where the magbar colour scale is placed within the plot. Specify one of 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'bottomright' and 'centre'.
}
  \item{\dots}{
Arguments to pass to magplot (hyper.plot2d) or plot3d (hyper.plot3d). When executing the hyper.fit class function plot.hyper.fit dots/ellipses are first passed to hyper.plot2d / hyper.plot3d and then onto magplot / plot3d for any unmatched arguments.
}
}
\value{
The plotting functions also return the sigma tension of the data points given the inputs. i.e. the output of hyper.like with output='sig'.
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
#### A very simple 2D example ####

#Make the simple data:

simpledata=cbind(x=1:10,y=c(12.2, 14.2, 15.9, 18.0, 20.1, 22.1, 23.9, 26.0, 27.9, 30.1))
simpfit=hyper.fit(simpledata)
summary(simpfit)
plot(simpfit)

#Increase the scatter:

simpledata2=cbind(x=1:10,y=c(11.6, 13.7, 15.5, 18.2, 21.2, 21.5, 23.6, 25.6, 27.9, 30.1))
simpfit2=hyper.fit(simpledata2)
summary(simpfit2)
plot(simpfit2)

#Assuming the error in each y data point is the same sy=0.5, we no longer need any
#component of intrinsic scatter to explain the data:

simpledata2err=cbind(sx=0, sy=rep(0.5, length(simpledata2[, 1])))
simpfit2werr=hyper.fit(simpledata2, vars=simpledata2err)
summary(simpfit2werr)
plot(simpfit2werr)

#### Simple Example in hyper.fit paper ####

#Fit with no error:

xval = c(-1.22, -0.78, 0.44, 1.01, 1.22)
yval = c(-0.15, 0.49, 1.17, 0.72, 1.22)

fitnoerror=hyper.fit(cbind(xval, yval))
plot(fitnoerror)

#Fit with independent x and y error:

xerr = c(0.12, 0.14, 0.20, 0.07, 0.06)
yerr = c(0.13, 0.03, 0.07, 0.11, 0.08)
fitwitherror=hyper.fit(cbind(xval, yval), vars=cbind(xerr, yerr)^2)
plot(fitwitherror)

#Fit with correlated x and y error:

xycor = c(0.90, -0.40, -0.25, 0.00, -0.20)
fitwitherrorandcor=hyper.fit(cbind(xval, yval), covarray=makecovarray2d(xerr, yerr, xycor))
plot(fitwitherrorandcor)

#### A 2D example with fitting a line ####

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
#Or even easier:
plot(fitline, trans=0.2, asp=1)

#### A 3D example with fitting a plane ####

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
fitplane=hyper.fit(X=X, covarray=covarray, coord.type='theta', scat.type='orth')
hyper.plot3d(X=X, covarray=covarray, fitobj=fitplane)
#Or even easier:
plot(fitplane)
}

#### Example using the data from Hogg 2010 ####

#Example using the data from Hogg 2010: http://arxiv.org/pdf/1008.4686v1.pdf

#Full data

\dontrun{

data(hogg)
fithogg=hyper.fit(X=cbind(hogg$x, hogg$y), covarray=makecovarray2d(hogg$x_err, hogg$y_err,
hogg$corxy), coord.type='theta', scat.type='orth')
hyper.plot2d(X=cbind(hogg$x, hogg$y), covarray=makecovarray2d(hogg$x_err, hogg$y_err,
hogg$corxy), fitobj=fithogg, trans=0.2, xlim=c(0, 300), ylim=c(0, 700))
#Or even easier:
plot(fithogg, trans=0.2)

#We now do exercise 17 of Hogg 2010 using trimmed data, where we remove the high tension
#data point 3 (which we can see as the reddest point in the above plot:

data(hogg)
hoggtrim=hogg[-3,]
fithoggtrim=hyper.fit(X=cbind(hoggtrim$x, hoggtrim$y), covarray=makecovarray2d(hoggtrim$x_err,
hoggtrim$y_err, hoggtrim$corxy), coord.type='theta', scat.type='orth', algo.func='LA')
hyper.plot2d(X=cbind(hoggtrim$x, hoggtrim$y), covarray=makecovarray2d(hoggtrim$x_err,
hoggtrim$y_err, hoggtrim$corxy), fitobj=fithoggtrim, trans=0.2, xlim=c(0, 300), ylim=c(0, 700))
#Or even easier:
plot(fithoggtrim, trans=0.2)

#We can compare this against the previous fit with:
hyper.plot2d(cbind(hoggtrim$x, hoggtrim$y), covarray=makecovarray2d(hoggtrim$x_err,
hoggtrim$y_err, hoggtrim$corxy), fitobj=fithogg, trans=0.2, xlim=c(0, 300), ylim=c(0, 700))

}

#### Example using 'real' data with intrinsic scatter ####

\dontrun{

data(intrin)
fitintrin=hyper.fit(X=cbind(intrin$x, intrin$y), vars=cbind(intrin$x_err,
intrin$y_err)^2, coord.type='theta', scat.type='orth', algo.func='LA')
hyper.plot2d(cbind(intrin$x, intrin$y), covarray=makecovarray2d(intrin$x_err,
intrin$y_err, intrin$corxy), fitobj=fitintrin, trans=0.1, pch='.', asp=1)
#Or even easier:
plot(fitintrin, trans=0.1, pch='.', asp=1)

}

#### Example using flaring trumpet data ####

\dontrun{

data(trumpet)
fittrumpet=hyper.fit(X=cbind(trumpet$x, trumpet$y), covarray=makecovarray2d(trumpet$x_err,
trumpet$y_err, trumpet$corxy), coord.type='theta', algo.func='LA')
hyper.plot2d(cbind(trumpet$x, trumpet$y), covarray=makecovarray2d(trumpet$x_err,
trumpet$y_err, trumpet$corxy), fitobj=fittrumpet, trans=0.1, pch='.', asp=1)
#Or even easier:
plot(fittrumpet, trans=0.1, pch='.', asp=1)
#If you look at the ?hyper.fit example we find that zero intrinsic scatter is actually
#preferred, but we don't see this in the above plot.

}

#### Example using 6dFGS Fundamental Plane data ####

\dontrun{

data(FP6dFGS)
fitFP6dFGSw=hyper.fit(FP6dFGS[,c('logIe_J', 'logsigma', 'logRe_J')],
vars=FP6dFGS[,c('logIe_J_err', 'logsigma_err', 'logRe_J_err')]^2, weights=FP6dFGS[,'weights'],
coord.type='alpha', scat.type='vert.axis')
#We turn the ellipse plotting off to speed things up:
plot(fitFP6dFGSw, doellipse=FALSE, alpha=0.5)

}

#### Example using GAMA mass-size relation data ####

\dontrun{

data(GAMAsmVsize)
fitGAMAsmVsize=hyper.fit(GAMAsmVsize[,c('logmstar', 'rekpc')],
vars=GAMAsmVsize[,c('logmstar_err', 'rekpc_err')]^2, weights=GAMAsmVsize[,'weights'],
coord.type='alpha', scat.type='vert.axis')
#We turn the ellipse plotting off to speed things up:
plot(fitGAMAsmVsize, doellipse=FALSE, unlog='x')

#This is obviously a poor fit since the y data has a non-linear dependence on x. Let's try
#using the logged y-axis and converted errors:

fitGAMAsmVsizelogre=hyper.fit(GAMAsmVsize[,c('logmstar', 'logrekpc')],
vars=GAMAsmVsize[,c('logmstar_err', 'logrekpc_err')]^2, weights=GAMAsmVsize[,'weights'],
coord.type='alpha', scat.type='vert.axis')
#We turn the ellipse plotting off to speed things up:
plot(fitGAMAsmVsizelogre, doellipse=FALSE, unlog='xy')

#We can compare to a fit with no errors used:

fitGAMAsmVsizelogrenoerr=hyper.fit(GAMAsmVsize[,c('logmstar', 'logrekpc')],
weights=GAMAsmVsize[,'weights'], coord.type='alpha', scat.type='vert.axis')
#We turn the ellipse plotting off to speed things up:
plot(fitGAMAsmVsizelogrenoerr, doellipse=FALSE, unlog='xy')

}

### Example using Tully-Fisher relation data ###

\dontrun{

data(TFR)
TFRfit=hyper.fit(X=TFR[,c('logv','M_K')],vars=TFR[,c('logv_err','M_K_err')]^2)
plot(TFRfit, xlim=c(1.7,2.5), ylim=c(-19,-26))

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ plot }
\keyword{ hyper }
\keyword{ linear }
\keyword{ plane }
\keyword{ regression }
