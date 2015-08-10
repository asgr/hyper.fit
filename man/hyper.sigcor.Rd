\name{hyper.sigcor}
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
  Robotham, A.S.G., & Obreschkow, D., PASA, in press
}
\author{
  Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.basic}}, \code{\link{hyper.convert}}, \code{\link{hyper.data}}, \code{\link{hyper.fit}}, \code{\link{hyper.plot}}, \code{\link{hyper.sigcor}}, \code{\link{hyper.summary}}
}
\examples{
#The below will take *a long* time to run- of the order a few days for the LD tests.
\dontrun{
Ngen=1e3
sdsamp=3
testvec=c(5,10,20,50,100)

set.seed(650)

fittestmean={}
for(Nsamp in testvec){
  print(paste('Nsamp=', Nsamp))
  fittest=matrix(0, nrow=Ngen, ncol=3)
  for(i in 1:Ngen){
    if(i %% 100==0){print(paste('Ngen=',i))}
    mockx=runif(Nsamp, -100, 100)
    mocky=mockx*tan(45*pi/180)+rnorm(Nsamp, sd=sdsamp)
    fittest[i,]=hyper.fit(X=cbind(mockx,mocky), coord.type='theta', scat.type='vert.axis')$parm
  }
  convtest2dOpt=rbind(convtest2dOpt, c(N=Nsamp,Raw=mean(fittest[,3]),
  mean(fittest[,3])*hyper.sigcor(Nsamp, 2)))
}

fittestmeanLD={}
for(Nsamp in testvec){
  print(paste('Nsamp=', Nsamp))
  fittest=matrix(0, nrow=Ngen, ncol=3)
  for(i in 1:Ngen){
    if(i %% 100==0){print(paste('Ngen=',i))}
    mockx=runif(Nsamp, -100, 100)
    mocky=mockx*tan(45*pi/180)+rnorm(Nsamp, sd=sdsamp)
    fittest[i,]=hyper.fit(X=cbind(mockx,mocky), coord.type='theta', scat.type='vert.axis',
    algo.func='LD', algo.method='GG', Specs=list(Grid=seq(-0.1,0.1, len=5), dparm=NULL,
    CPUs=1, Packages=NULL, Dyn.libs=NULL))$parm
    print(fittest[i,])
  }
  convtest2dLD=rbind(convtest2dLD, c(N=Nsamp, Raw=mean(fittest[,3]),
  mean(fittest[,3])*hyper.sigcor(Nsamp, 2)))
}

normtestmean={}
for(Nsamp in testvec){
  print(paste('Nsamp=', Nsamp))
  normtest={}
  for(i in 1:Ngen){
    if(i %% 100==0){print(paste('Ngen=',i))}
    normtemp=rnorm(Nsamp, sd=sdsamp)
    normtest=c(normtest, sqrt(sum((normtemp-mean(normtemp))^2)/Nsamp))
  }
  convtest1dNorm=rbind(convtest1dNorm, c(N=Nsamp, Raw=mean(normtest),
  mean(normtest)*hyper.sigcor(Nsamp, 1)))
}
}
#The runs above have been pre-generated and can be loaded via

data(convtest2dOpt)
data(convtest2dLD)
data(convtest1dNorm)

magplot(convtest2dOpt[,c('N','Raw')],xlim=c(5,100),ylim=c(0,4),type='b',log='x')
lines(convtest2dOpt[,c('N','sampbias2popunbias')],type='b',lty=2,pch=4)
lines(convtest2dLD[,c('N','Raw')],type='b',col='blue')
lines(convtest2dLD[,c('N','bias2unbias')],type='b',lty=2,pch=4,col='blue')
lines(convtest1dNorm[,c('N','Raw')],type='b',col='red')
lines(convtest1dNorm[,c('N','sampbias2popunbias')],type='b',lty=2,pch=4,col='red')
legend('topleft', legend=c('2 DoF and optim fit','2 DoF and LD fit', '1 DoF and direct SD'),
col=c('black','blue','red'),pch=1)
legend('topright', legend=c('Raw intrinsic scatter', 'Corrected intrinsic scatter'),
lty=c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ utility }
\keyword{ sigma }
\keyword{ bias }
\keyword{ sample }
\keyword{ population }
