\name{Data}
\alias{hogg}
\alias{intrin}
\alias{trumpet}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data included in hyper.fit package
}
\description{
hogg:     Toy ASCII table taken from Hogg 2010. Contains x/y/sx/sy/corxy columns. Provided by David Hogg: http://arxiv.org/pdf/1008.4686v1.pdf.

intrin:   Toy ASCII table that has intrinsic scatter. Contains x/y/sx/sy/corxy columns. Provided by Michelle Cluver.

trumpet:  Toy ASCII table that has trumpet-like covariance errors and no intrinsic scatter. Contains x/y/sx/sy/corxy columns. Provided by Johannes Buchner.
}
\usage{
data(hogg)
data(intrin)
data(trumpet)
}
\references{
Robotham, A.S.G., & Obreschkow, D., 2014

Hogg, D., Bovy, J., Lang, D., 2010 (http://arxiv.org/pdf/1008.4686v1.pdf)
}
\author{
Aaron Robotham and Danail Obreschkow
}
\seealso{
\code{\link{hyper.fit}}, \code{\link{hyper.plot2d}}, \code{\link{hyper.plot3d}}, \code{\link{makecovmat2d}}, \code{\link{makecovmat3d}}
}
\examples{
hogg=read.table(system.file('data/hogg.tab', package='hyper.fit'),header=TRUE)
print(hogg[1:10,])

intrin=read.table(system.file('data/intrin.tab', package='hyper.fit'), header=TRUE)
print(intrin[1:10,])

trumpet=read.table(system.file('data/trumpet.tab', package='hyper.fit'), header=TRUE)
print(trumpet[1:10,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ data }
