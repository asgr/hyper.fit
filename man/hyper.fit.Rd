\name{hyper.fit}
\alias{hyper.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Top level function that attempts to fit a hyperplane to provided data.
}
\description{
Top level fitting function that uses downhill searches (optim/LaplaceApproximation) or MCMC (LaplacesDemon) to search out the best fitting parameters for a hyperplane (minimum a 1D line for 2D data), including the intrinsic scatter as part of the fit.
}
\usage{
hyper.fit(X, covarray, vars, parm, parm.coord, parm.beta, parm.scat, parm.errorscale=1,
vert.axis, weights, k.vec, itermax = 1e4, coord.type = 'alpha', scat.type = 'vert.axis',
algo.func = 'optim', algo.method = 'default', Specs = list(Grid=seq(-0.1,0.1, len=5),
dparm=NULL, CPUs=1, Packages=NULL, Dyn.libs=NULL), doerrorscale = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A position matrix with N (number of data points) rows by d (number of dimensions) columns.
}
  \item{covarray}{
A dxdxN array containing the full covariance (d=dimensions, N=number of dxd matrices in the array stack). The makecovarray2d and makecovarray3d are convenience functions that make populating 2x2xN and 3x3xN arrays easier for a novice user.
}
  \item{vars}{
A variance matrix with the N (numver of data points) rows by d (number of dimensions) columns. In effect this is the diagonal elements of covarray where all other terms are zero. If covarray is also provided that input argument is used instead.
}
  \item{parm}{
Vector of all initial parameters. This should be a concatenation of c(parm.coord, parm.beta, parm.scat). If doerrorscale=TRUE then there should be an extrat element to give c(parm.coord, parm.beta, parm.scat, parm.errorscale). See the parm.coord, parm.beta, parm.scat and parm.errorscale arguments below for more information on what these different quantities describe.
}
  \item{parm.coord}{
Vector of initial coord parameters. These are either angles that produce the vectors that predict the vert.axis dimension (coord.type="theta"), the gradients of these (coord.type="alpha") or they are the vector elements normal to the hyperplane (coord.type="normvec"). Depending on the argument used, either beta along the vertical axis is generated (coord.type="alpha"/theta") or no beta is required (coord.type="normvec").
}
  \item{parm.beta}{
Initial value of beta. This is either specified as the absolute distance from the origin to the hyperplane or as the intersection of the hyperplane on the vert.axis dimension being predicted.
}
  \item{parm.scat}{
Initial value of the intrinsic scatter. This is either specified as the scatter orthogonal to the hyperplane or as the scatter along the vert.axis dimension being predicted.
}
  \item{parm.errorscale}{
If doerrorscale=TRUE, this argument is the initial errorscale (default is 1). See the doerrorscale argument below for more details.
}
  \item{vert.axis}{
Which axis should the hyperplane equation be formulated for. This must be a number which specifies the column of position matrix X to be defined by the hyperplane. If missing, then the projection dimension is assumed to be the last column of the X matrix.
}
  \item{weights}{
Vector of multiplicative weights for each row of the X data matrix. i.e. if this is 2 then it is equivalent to having two identical data points with weights equal to 1. Should be either of length 1 (in which case elements are repeated as required) or the same length as the number of rows in the data matrix X.
}
\item{k.vec}{
A vector defining the direction of an exponential sampling distribution in the data. The length is the scaling "a" of the generative exponent (i.e., exp(a*x)), and it points in the direction of *increasing* density (see example below). If provided, k.vec must be the same length as the dimensions of the data. k.vec has the most noticeable effect on the beta offset parameter. It is correcting for the phenomenom Astronomers often call Eddington bias. For discussion on the use of k.vec see Appendix B of Robotham & Obreschkow.
}
  \item{itermax}{
The maximum iterations to use for either the LaplaceApproximation function or LaplacesDemon function.
}
  \item{coord.type}{
This specifies whether the fit should be done in terms of the normal vector to the hyperplane (coord.type="normvec") gradients defined to produce values along the vert.axis dimension (coord.type="alpha") or by the values of the angles that form the gradients (coord.type="theta"). "alpha" is the default since it is the most common means of parameterising hyperplane fits in the astronomy literature.
}
  \item{scat.type}{
This specifies whether the intrinsic scatter should be defined orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis).
}
  \item{algo.func}{
If \option{algo.func="optim"} (default) hyper.fit will optimise using the R base \code{\link{optim}} function. If \option{algo.func="LA"} will optimise using the \code{\link{LaplaceApproximation}} function. If \option{algo.func="LD"} will optimise using the \code{\link{LaplacesDemon}} function. For both \option{algo.func="LA"} and \option{algo.func="LD"} the \code{LaplacesDemon} package will be used (see http://www.bayesian-inference.com/software).
}
  \item{algo.method}{
Specifies the \option{method} argument of \code{\link{optim}} function when using \option{algo.func="optim"} (if not specified hyper.fit will use "Nelder-Mead" for \code{optim}). Specifies \option{Method} argument of \code{\link{LaplaceApproximation}} function when using \option{algo.func="LA"} (if not specified hyper.fit will use "NM" for \code{LaplaceApproximation}). Specifies \option{Algorithm} argument of \code{\link{LaplacesDemon}} function when using \option{algo.func="LD"} (if not specified hyper.fit will use "CHARM" for \code{LaplacesDemon}). When using \option{algo.func="LD"} the user can also specify further options via the Specs argument below.
}
  \item{Specs}{
Inputs to pass to the LaplacesDemon function. Default Specs=list(alpha.star = 0.44) option is for the default CHARM algorithm (see algo.method above).
}
  \item{doerrorscale}{
If FALSE then the provided covariance are treated as it. If TRUE then the likelihood function is also allowed to rescale all errors by a uniform multiplicative value.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{optim}}, \code{\link{LaplaceApproximation}} or \code{\link{LaplacesDemon}} depending on the option used for \option{algo.func}.
}
}
\details{
Setting doerrorscale to TRUE allows for stable solutions when errors are overestimated, e.g. when the intrinsic scatter is equal to zero but the data is more clustered around the optimal likelihood plane than expected from the data covariance array. See Examples below for a 2D scenario where this is helpful.

algo.func="LD" also returns the probability that the generative model has exactly zero intrinsic scatter (zeroscatprob in the output). The other available functions (algo.func="optim" and algo.func="LA") can find exact solutions equal to zero since they are strictly mode finding (i.e. maximum likelihood) routines. algo.func="LD" is MCMC based, so naturally returns a mean/expectation which *must* have a finite positive value for the intrinsic scatter (it's not allowed to travel below zero). The zeroscatprob output works better with a Metropolis type scheme, e.g. algo.method='CHARM', Specs=list(alpha.star = 0.44) works well for many cases.
}
\value{

The function returns a multi-component list containing:

\item{parm}{
  Vector of the main paramter fit outputs specified as set by the coord.type and scat.type options.
}
\item{parm.vert.axis}{
  Vector of the main parameter fit outputs specified strictly along the last column dimension of X (for both the intrinsic scatter and the beta offset). The order of the columns in X therefore affects the contents of this vector. This output assume a coord.type="alpha" and scat.type="vert.axis".
}
\item{fit}{
  The direct output of the specified algo.func. So either the natural return from optim (class type "optim"), LaplaceApproximation (class type "laplace") or LaplacesDemon (class type "demonoid").
}
\item{sigcor}{
  If algo.func="optim" or algo.func="LA" then this list element contains the unbiased population estimator for the intrinsic scatter corrected via the sampbias2popunbias output of \code{\link{hyper.sigcor}}, i.e. it uses the full correction from hyper.sigcor (element 3 of the correction vector generated). If algo.func="LD" then this list element contains the unbiased population estimator for the intrinsic scatter corrected via the bias2unbias output of \code{\link{hyper.sigcor}}, i.e. it uses the standard deviation bias correction from hyper.sigcor (element 1 of the correction vector generated). See \code{\link{hyper.sigcor}} for details on the different outputs and the convergence tests the demonstrate why these different correction methods are desirable.
}
\item{parm.covar}{
  The covariance matrix for parm. Only for algo.func="optim" and algo.func="LA".
}
\item{zeroscatprob}{
  The fraction of samples for the intrinsic scatter which are at *exactly* zero, which provides a guideline probability for the intrinsic scatter being truly zero, rather than the expectation which will always be a finite amount above zero. Only for algo.func="LD". See Details above.
}
\item{hyper.plane}{
  Object class of type hyper.plane.param, as output by hyper.covert. This standardises the final fit in the users requested coord.type and scat.type, and allows easy conversion to other systems by using the class dependent convert function on it.
}
\item{N}{
  The number of rows of matrix X.
}
\item{dims}{
  The number of columns of matrix X.
}
\item{X}{
  The input matrix X.
}
\item{covarray}{
  The covarray used for fitting.
}
\item{weights}{
  The weights used for fitting.
}
\item{call}{
  The actual call used to run hyper.fit.
}
\item{args}{
  The arguments used to run hyper.fit.
}
\item{LL}{
  A list containing elements sum (the total log-likelihodd), val (the individual log likelihoods) and sig (the effective sigma offsets). See hyper.like for more information on these outputs.
}
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

#We can fit for 6 different combinations of coordinate system:

print(hyper.fit(simpledata, coord.type='theta', scat.type='orth')$parm)
print(hyper.fit(simpledata, coord.type='alpha', scat.type='orth')$parm)
print(hyper.fit(simpledata, coord.type='normvec', scat.type='orth')$parm)
print(hyper.fit(simpledata, coord.type='theta', scat.type='vert.axis')$parm)
print(hyper.fit(simpledata, coord.type='alpha', scat.type='vert.axis')$parm)
print(hyper.fit(simpledata, coord.type='normvec', scat.type='vert.axis')$parm)

#These all describe the same hyperplane (or line in this case). We can convert between
#systems by using the hyper.convert utility function:

fit4normvert=hyper.fit(simpledata, coord.type='normvec', scat.type='vert.axis')$parm
hyper.convert(fit4normvert, in.coord.type='normvec', out.coord.type='theta',
in.scat.type='vert.axis', out.scat.type='orth')$parm

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

\dontrun{

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
plot(fitline, trans=0.2, asp=1)

#We can test to see if the errors are compatable with the intrinsic scatter:

fitlineerrscale=hyper.fit(X=X, covarray=covarray, coord.type='theta', doerrorscale=TRUE)
summary(fitlineerrscale)
plot(fitline, parm.errorscale=fitlineerrscale$parm['errorscale'], trans=0.2, asp=1)

#Within errors the errorscale parameter is 1, i.e. the errors are realistic, which we know
#they should be a priori since we made them ourselves.

}

#### A 2D example with exponential sampling & fitting a line ####

#Setup the initial data:

\dontrun{

set.seed(650)

#The effect of an exponential density function along y is to offset the Gaussian mean by
#0.5 times the factor 'a' in exp(a*x), i.e.:

normfac=dnorm(0,sd=1.1)/(dnorm(10*1.1^2,sd=1.1)*exp(10*10*1.1^2))
magplot(seq(5,15,by=0.01), normfac*dnorm(seq(5,15, by=0.01), sd=1.1)*exp(10*seq(5,15,
by=0.01)), type='l')
abline(v=10*1.1^2,lty=2)

#The above will not be correctly normalised to form a true PDF, but the shift in the mean
#is clear, and it doesn't alter the standard deviation at all:

points(seq(5,15,by=0.1), dnorm(seq(5,15, by=0.1), mean=10*1.1^2, sd=1.1),col='red')

#Applying the same principal to our random data we apply the offset due to our exponential
#generative slope in y:

set.seed(650)

sampN=200
vert.scat=10
sampexp=0.1
ang=30
randatax=runif(200,-100,100)
randatay=randatax*tan(ang*pi/180)+rnorm(sampN, mean=sampexp*vert.scat^2, sd=vert.scat)
sx=runif(sampN, 0, 10); sy=runif(sampN, 0, 10)

mockvararray=makecovarray2d(sx, sy, corxy=0)

errxy={}
for(i in 1:sampN){
  rancovmat=ranrotcovmat2d(mockvararray[,,i])
  errxy=rbind(errxy, mvrnorm(1, mu=c(0, sampexp*sy[i]^2), Sigma=rancovmat))
  mockvararray[,,i]=rancovmat
  }
randatax=randatax+errxy[,1]
randatay=randatay+errxy[,2]
sx=sqrt(mockvararray[1,1,]); sy=sqrt(mockvararray[2,2,]); corxy=mockvararray[1,2,]/(sx*sy)
mock=data.frame(x=randatax, y=randatay, sx=sx, sy=sy, corxy=corxy)

#Do the fit. Notice that the second element of k.vec has the positive sign, i.e. we are moving
#data that has been shifted positively by the positive exponential slope in y back to where it
#would exist without the slope (i.e. if it had an equal chance of being scattered in both
#directions, rather than being preferentially offset in the direction away from denser data).
#This dense -> less-dense shift i s known as Eddington bias in astronomy, and is common in all
#power-law distributions that have intrinsic scatter (e.g. Schechter LF and dark matter HMF).

X=cbind(mock$x, mock$y)
covarray=makecovarray2d(mock$sx, mock$sy, mock$corxy)
fitlineexp=hyper.fit(X=X, covarray=covarray, coord.type='theta', k.vec=c(0,sampexp),
scat.type='vert.axis')
summary(fitlineexp)
plot(fitlineexp, k.vec=c(0,sampexp))

#If we ignore the k.vec when calculating the plotting sigma values you can see it has
#a significant effect:

plot(fitlineexp, trans=0.2, asp=1)

#Compare this to not including the known exponential slope:

fitlinenoexp=hyper.fit(X=X, covarray=covarray, coord.type='theta', k.vec=c(0,0),
scat.type='vert.axis')
summary(fitlinenoexp)
plot(fitlinenoexp, trans=0.2, asp=1)

}

#The theta and intrinsic scatter are similar, but the offset is shifted significantly
#away from zero.

#### A 3D example with fitting a plane ####

#Setup the initial data:

\dontrun{

set.seed(650)
sampN=200
initscat=3
randatax=runif(sampN, -100, 100)
randatay=runif(sampN, -100, 100)
randataz=rnorm(sampN, sd=initscat)
sx=runif(sampN, 0, 5); sy=runif(sampN, 0, 5); sz=runif(sampN, 0, 5)

mockvararray=makecovarray3d(sx, sy, sz, corxy=0, corxz=0, coryz=0)

errxyz={}
for(i in 1:sampN){
  rancovmat=ranrotcovmat3d(mockvararray[,,i])
  errxyz=rbind(errxyz,mvrnorm(1, mu=c(0, 0, 0), Sigma=rancovmat))
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

X=cbind(mockplane$x, mockplane$y, mockplane$z)
covarray=makecovarray3d(mockplane$sx, mockplane$sy, mockplane$sz, mockplane$corxy,
mockplane$corxz, mockplane$coryz)
fitplane=hyper.fit(X=X, covarray=covarray, coord.type='theta', scat.type='orth')
summary(fitplane)
plot(fitplane)

}

#### Example using the data from Hogg 2010 ####

#Example using the data from Hogg 2010: http://arxiv.org/pdf/1008.4686v1.pdf

#Load data:

\dontrun{

data(hogg)

#Fit:

fithogg=hyper.fit(X=cbind(hogg$x, hogg$y), covarray=makecovarray2d(hogg$x_err, hogg$y_err,
hogg$corxy), coord.type='theta', scat.type='orth')
summary(fithogg)
plot(fithogg, trans=0.2)

#We now do exercise 17 of Hogg 2010 using trimmed data:

hoggtrim=hogg[-3,]
fithoggtrim=hyper.fit(X=cbind(hoggtrim$x, hoggtrim$y), covarray=makecovarray2d(hoggtrim$x_err,
hoggtrim$y_err, hoggtrim$corxy), coord.type='theta', scat.type='orth', algo.func='LA')
summary(fithoggtrim)
plot(fithoggtrim, trans=0.2)

#We can get more info from looking at the Summary1 output of the LaplaceApproximation:

print(fithoggtrim$fit$Summary1)

#MCMC (exercise 18):

fithoggtrimMCMC=hyper.fit(X=cbind(hoggtrim$x, hoggtrim$y), covarray=
makecovarray2d(hoggtrim$x_err, hoggtrim$y_err, hoggtrim$corxy), coord.type='theta',
scat.type='orth', algo.func='LD', algo.method='CHARM', Specs=list(alpha.star = 0.44))
summary(fithoggtrimMCMC)

#We can get additional info from looking at the Summary1 output of the LaplacesDemon:

print(fithoggtrimMCMC$fit$Summary2)

magplot(density(fithoggtrimMCMC$fit$Posterior2[,3]), xlab='Intrinsic Scatter',
ylab='Probability Density')
abline(v=quantile(fithoggtrimMCMC$fit$Posterior2[,3], c(0.95,0.99)), lty=2)

}

#### Example using 'real' data with intrinsic scatter ####

\dontrun{

data(intrin)

fitintrin=hyper.fit(X=cbind(intrin$x, intrin$y), vars=cbind(intrin$x_err,
intrin$y_err)^2, coord.type='theta', scat.type='orth', algo.func='LA')
summary(fitintrin)
plot(fitintrin, trans=0.1, pch='.', asp=1)

fitintrincor=hyper.fit(X=cbind(intrin$x, intrin$y), covarray=makecovarray2d(intrin$x_err,
intrin$y_err, intrin$corxy), coord.type='theta', scat.type='orth', algo.func='LA')
summary(fitintrincor)
plot(fitintrincor, trans=0.1, pch='.', asp=1)

}

#### Example using flaring trumpet data ####

\dontrun{

data(trumpet)
fittrumpet=hyper.fit(X=cbind(trumpet$x, trumpet$y), covarray=makecovarray2d(trumpet$x_err,
trumpet$y_err, trumpet$corxy), coord.type='normvec', algo.func='LA')
summary(fittrumpet)
plot(fittrumpet, trans=0.1, pch='.', asp=1)

#The best fit solution has a scat.orth very close to 0, so it is worth considering if the
#data should truly have 0 intrinsic scatter.

}

\dontrun{

#To find the likelihood of zero intrinsic scatter we will need to run LaplacesDemon. The
#following will take a couple of minutes to run:

set.seed(650)
fittrumpetMCMC=hyper.fit(X=cbind(trumpet$x, trumpet$y), covarray=makecovarray2d(trumpet$x_err,
trumpet$y_err, trumpet$corxy), coord.type='normvec', algo.func='LD', itermax=1e5)

#Assuming the user has specified the same initial seed we should find that the data
#has exactly zero intrinsic scatter with ~47\% likelihood:

print(fittrumpetMCMC$zeroscatprob)

#We can also make an assessment of whether the data has even less scatter than expected
#given the expected errors:

set.seed(650)
fittrumpetMCMCerrscale=hyper.fit(X=cbind(trumpet$x, trumpet$y), covarray=makecovarray2d(
trumpet$x_err, trumpet$y_err, trumpet$corxy), itermax=1e5, coord.type='normvec', algo.func='LD',
algo.method='CHARM', Specs=list(alpha.star = 0.44), doerrorscale=TRUE)

#Assuming the user has specified the same initial seed we should find that the data
#has exactly zero intrinsic scatter with ~69\% likelihood:

print(fittrumpetMCMCerrscale$zeroscatprob)

}

#### Example using 6dFGS Fundamental Plane data ####

\dontrun{

data(FP6dFGS)

#First we try the fit without using any weights:

fitFP6dFGS=hyper.fit(FP6dFGS[,c('logIe_J', 'logsigma', 'logRe_J')], 
vars=FP6dFGS[,c('logIe_J_err', 'logsigma_err', 'logRe_J_err')]^2, coord.type='alpha',
scat.type='vert.axis')
summary(fitFP6dFGS)
plot(fitFP6dFGS, doellipse=FALSE, alpha=0.5)

}

#Next we add the censoring weights provided by C. Magoulas:

\dontrun{

fitFP6dFGSw=hyper.fit(FP6dFGS[,c('logIe_J', 'logsigma', 'logRe_J')],
vars=FP6dFGS[,c('logIe_J_err', 'logsigma_err', 'logRe_J_err')]^2, weights=FP6dFGS[,'weights'],
coord.type='alpha', scat.type='vert.axis')
summary(fitFP6dFGSw)
plot(fitFP6dFGSw, doellipse=FALSE, alpha=0.5)

#It is interesting to note the scatter orthogonal to the plane for the fundmental plane:

print(hyper.convert(coord=fitFP6dFGSw$parm[1:2], beta=fitFP6dFGSw$parm[3],
scat=fitFP6dFGSw$parm[4], in.scat.type='vert.axis', out.scat.type='orth',
in.coord.type='alpha'))

}

#### Example using GAMA mass-size relation data ####

\dontrun{

data(GAMAsmVsize)
fitGAMAsmVsize=hyper.fit(GAMAsmVsize[,c('logmstar', 'rekpc')],
vars=GAMAsmVsize[,c('logmstar_err', 'rekpc_err')]^2, weights=GAMAsmVsize[,'weights'],
coord.type='alpha', scat.type='vert.axis')
summary(fitGAMAsmVsize)
#We turn the ellipse plotting off to speed things up:
plot(fitGAMAsmVsize, doellipse=FALSE, unlog='x')

#This is obviously a poor fit since the y data has a non-linear dependence on x. Let's try
#using the logged y-axis and converted errors:

fitGAMAsmVsizelogre=hyper.fit(GAMAsmVsize[,c('logmstar', 'logrekpc')],
vars=GAMAsmVsize[,c('logmstar_err', 'logrekpc_err')]^2, weights=GAMAsmVsize[,'weights'],
coord.type='alpha', scat.type='vert.axis')
summary(fitGAMAsmVsizelogre)
#We turn the ellipse plotting off to speed things up:
plot(fitGAMAsmVsizelogre, doellipse=FALSE, unlog='xy')

#We can compare to a fit with no errors used:

fitGAMAsmVsizelogrenoerr=hyper.fit(GAMAsmVsize[,c('logmstar', 'logrekpc')],
weights=GAMAsmVsize[,'weights'], coord.type='alpha', scat.type='vert.axis')
summary(fitGAMAsmVsizelogrenoerr)
#We turn the ellipse plotting off to speed things up:
plot(fitGAMAsmVsizelogrenoerr, doellipse=FALSE, unlog='xy')

}

### Example using Tully-Fisher relation data ###

\dontrun{

data(TFR)
TFRfit=hyper.fit(X=TFR[,c('logv','M_K')],vars=TFR[,c('logv_err','M_K_err')]^2)
plot(TFRfit, xlim=c(1.7,2.5), ylim=c(-19,-26))

}

### Mase-Angular Momentum-Bulge/Total ###

\dontrun{

data(MJB)
MJBfit=hyper.fit(X=MJB[,c('logM','logj','B.T')], covarray=makecovarray3d(MJB$logM_err,
MJB$logj_err, MJB$B.T_err, MJB$corMJ, 0, 0))
plot(MJBfit)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ hyper }
\keyword{ linear }
\keyword{ plane }
\keyword{ regression }
