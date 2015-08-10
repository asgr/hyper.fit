\name{hyper.data}
\alias{hyper.data}
\alias{hogg}
\alias{intrin}
\alias{trumpet}
\alias{FP6dFGS}
\alias{GAMAsmVsize}
\alias{convtest2dOpt}
\alias{convtest2dLD}
\alias{convtest1dNorm}
\alias{TFR}
\alias{MJB}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data included in hyper.fit package
}
\description{
hogg:     Toy ASCII table taken from Hogg 2010. Contains x/y/sx/sy/corxy columns. Provided by David Hogg.

intrin:   Toy ASCII table that has intrinsic scatter. Contains x/y/sx/sy/corxy columns. Provided by Michelle Cluver.

trumpet:  Toy ASCII table that has trumpet-like covariance errors and no intrinsic scatter. Contains x/y/sx/sy/corxy columns. Provided by Johannes Buchner.

FP6dFGS: 6dFGS fundamental plane data taken from table 8 of Campbell et al 2014. We use columns 6/7, 12/13, 18/19 for the FP parameters and their errors, and column 26 = 111111 to create a clean selection. Provided By Christina Magoulas.

GAMAsmVsize: Galaxy mass vs size data taken from Lange et al 2014. Bottom-right panel of Figure 3 (i.e. r-band elliptical relation data). Provided by Rebecca Lange.

TFR: Tully-Fisher Relation data taken from Obreschkow and Meyer 2013. Provided by Danail Obreschkow.

MJB: Mase-Angular Momentum-Bulge/Total data taken from Obreschkow & Glazebrook 2014. Provided by Danail Obreschkow.

convtest2dOpt: Intrinsic scatter convergence test data for 2 DoF and optim fitted simulations, shown as an example in \code{\link{hyper.sigcor}}.

convtest2dLD: Intrinsic scatter convergence test data for 2 DoF and LaplacesDemon fitted simulations, shown as an example in \code{\link{hyper.sigcor}}.

convtest1dNorm: Intrinsic scatter convergence test data for 1 DoF and direct sample SD estimation, shown as an example in \code{\link{hyper.sigcor}}.
}
\usage{
data(hogg)
data(intrin)
data(trumpet)
data(FP6dFGS)
data(GAMAsmVsize)
data(TFR)
data(MJB)
data(convtest2dOpt)
data(convtest2dLD)
data(convtest1dNorm)
}
\references{
Robotham, A.S.G., & Obreschkow, D., PASA, in press

Campbell, L., et al., 2014, MNRAS, 443, 1231 (http://mnras.oxfordjournals.org/content/443/2/1231)

Cluver, M., et al., 2014, ApJ, 782, 90 (http://arxiv.org/pdf/1401.0837v1.pdf)

Hogg, D., Bovy, J., Lang, D., 2010 (http://arxiv.org/pdf/1008.4686v1.pdf)

Lange, R., et al., MNRAS, accepted

Obreschkow & Meyer, 2013, ApJ, 777, 140

Obreschkow & Glazebrook, 2014, ApJ, 784, 26
}
\author{
Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.basic}}, \code{\link{hyper.convert}}, \code{\link{hyper.data}}, \code{\link{hyper.fit}}, \code{\link{hyper.plot}}, \code{\link{hyper.sigcor}}, \code{\link{hyper.summary}}
}
\examples{
hogg=read.table(system.file('data/hogg.tab', package='hyper.fit'),header=TRUE)
#or
data(hogg)
print(hogg[1:10,])

intrin=read.table(system.file('data/intrin.tab', package='hyper.fit'), header=TRUE)
#or
data(intrin)
print(intrin[1:10,])

trumpet=read.table(system.file('data/trumpet.tab', package='hyper.fit'), header=TRUE)
#or
data(trumpet)
print(trumpet[1:10,])

FP6dFGS=read.table(system.file('data/FP6dFGS.tab', package='hyper.fit'), header=TRUE)
#or
data(FP6dFGS)
print(FP6dFGS[1:10,])

GAMAsmVsize=read.table(system.file('data/GAMAsmVsize.tab', package='hyper.fit'), header=TRUE)
#or
data(GAMAsmVsize)
print(GAMAsmVsize[1:10,])

TFR=read.table(system.file('data/TFR.tab', package='hyper.fit'), header=TRUE)
#or
data(TFR)
print(TFR[1:10,])

MBJ=read.table(system.file('data/MJB.tab', package='hyper.fit'), header=TRUE)
#or
data(MJB)
print(MJB[1:10,])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ data }
