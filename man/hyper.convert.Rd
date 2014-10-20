\name{hyper.convert}
\alias{hyper.convert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameterisation conversion functions.
}
\description{
This function allows the user to generically convert their current plane definition system to an alternative. The obvious use case might be when the user has an equation defined as a projection formula along a preferred axis (e.g. z=ax+by+c) and they want to find the orthogonal offset of the hyperplane to the origin, and the intrinsic scatter orthogonal to the hyperplane.
}
\usage{
hyper.convert(parm, coord, beta = 0, scat = 0, in.coord.type = "alpha", out.coord.type,
in.proj.type = "vert.axis", out.proj.type, in.vert.axis, out.vert.axis)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Vector of all parameters. This should be a concatenation of c(coord,beta,scat). Either 'parm' or 'coord' must be specified. 
}
  \item{coord}{
The current coordinate parameters using the in.coord.type coordinate system. This should be a vector of length dimensions-1 (i.e. 1 for 2D xy data and 2 for 3D xyz data). The coord argument must be explictly specified by the user.
}
  \item{beta}{
The current offset of the hyperplane defined using the in.proj.type projection system. The default is 0.
}
  \item{scat}{
The current intrinsic scatter of the hyperplane defined using the in.proj.type projection system. The default is 0.
}
  \item{in.coord.type}{
This specifies whether the input coord vector is defined in terms of the normal vector to the hyperplane (normvec) gradients defined to produce values along the vert.axis dimension (alpha, default) or by the values of the angles that form the gradients (theta).
}
  \item{out.coord.type}{
This specifies whether the output coord vector is defined in terms of the normal vector to the hyperplane (normvec) gradients defined to produce values along the vert.axis dimension (alpha) or by the values of the angles that form the gradients (theta). If missing it takes the value of in.coord.type.
}
  \item{in.proj.type}{
This specifies whether the input beta/scat are defined as orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis).
}
  \item{out.proj.type}{
This specifies whether the output beta/scat should be defined as orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis). If missing it takes the value of in.proj.type.
}
  \item{in.vert.axis}{
This specifies the input vertical projection axis. If missing it uses the maximum dimension value (i.e. y-axis for a 2d dataset).
}
  \item{out.vert.axis}{
This specifies the output/requested vertical projection axis. If missing it takes the value of in.vert.axis.
}
}

\value{

The function returns a multi-component list containing:
  \item{parm}{
  parm is a concatination of the main parameters, e.g. parm=c(coord, beta, scat). It can be used as a direct input to other hyper.fit functions that can accept a parm type input.
}
  \item{coord}{
  The output coordinate parameters using the out.coord.type coordinate system.
}
  \item{beta}{
  The output offset of the hyperplane defined using the out.proj.type projection system. Potentially standardised by the abs.beta.orth argument.
}
  \item{scat}{
  The output intrinsic scatter of the hyperplane defined using the out.proj.type projection system.
}
  \item{unitvec}{
  The unit vector orthogonal to the hyperplane. This is always the vector pointing from the origin *to* the hyperplane. 
}
  \item{beta.orth}{
  The absolute distance of the hyperplane to the origin. By definition this will always be positive.
}
  \item{scat.orth}{
  The intrinsic scatter orthogonal to the hyperplane. By definition this will always be positive.
}
  \item{coord.type}{
  The requested out.coord.type.
}
  \item{proj.type}{
  The requested out.proj.type.
}
  \item{vert.axis}{
  The requested out.vert.axis (if out.vert.axis is not specified, this will be the maximum dimension value, see in.vert.axis in arguments section above).  
}
}
\references{
Robotham, A.S.G., & Obreschkow, D., 2014
}
\author{
Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.basic}}, \code{\link{hyper.convert}}, \code{\link{hyper.data}}, \code{\link{hyper.fit}}, \code{\link{hyper.plot}}, \code{\link{hyper.sigcor}}, \code{\link{hyper.summary}}
}
\examples{

#Here we are assuming our plane forula is z=2x+3y+1 plus Gaussian intrinsic scatter along
#z with sd=4:

excoord.alpha=c(2,3)
exbeta.vert.axis=1
exscat.vert.axis=4

print(hyper.convert(coord=excoord.alpha, beta=exbeta.vert.axis, scat=exscat.vert.axis,
out.coord.type='theta', out.proj.type='orth'))

print(hyper.convert(coord=excoord.alpha, beta=exbeta.vert.axis, scat=exscat.vert.axis,
out.vert.axis=2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ convert }
\keyword{ hyper }% __ONLY ONE__ keyword per line
