\name{hyper.convert}
\alias{coord.convert}
\alias{beta.convert}
\alias{scat.convert}
\alias{vert.axis.convert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameterisation conversion functions.
}
\description{
These functions allow the user to generically convert their current plane definition systems to an alternative. The obvious use case might be when the user has an equation defined as projection formula along a preferred axis (e.g. z=ax+by+c) and they want to find the orthogonal offset of the hyperplane to the origin, and the intrinsic scatter orthogonal to the hyperplane.
}
\usage{
coord.convert(coord, in.coord.type = "alpha", out.coord.type = "alpha")

beta.convert(beta, coord, in.proj.type='vert.axis', out.proj.type='vert.axis',
in.coord.type='alpha')

scat.convert(scat,coord, in.proj.type='vert.axis', out.proj.type='vert.axis',
in.coord.type='alpha')

vert.axis.convert(coord, beta=0, scat=0, in.vert.axis, out.vert.axis,
in.proj.type='vert.axis', in.coord.type='alpha')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coord}{
The current coordinate parameters using the in.coord.type coordinate system. This should be a vector of length dimensions-1 (i.e. 1 for 2D xy data and 2 for 3D xyz data).
}
  \item{beta}{
The current offset of the hyper plane defined using the in.proj.type projection system.
}
  \item{scat}{
The current intrinsic scatter of the hyper plane defined using the in.proj.type projection system.
}
  \item{in.coord.type}{
This specifies whether the input coord parameter is defined in terms of the unit vector of the line (alpha) or by the values of the angles that form the unit vector (theta).
}
  \item{out.coord.type}{
This specifies whether the output coord parameter should be defined in terms of the unit vector of the line (alpha) or by the values of the angles that form the unit vector (theta).
}
  \item{in.proj.type}{
This specifies whether the input beta/scat are defined as orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis).
}
  \item{out.proj.type}{
This specifies whether the output beta/scat should be defined as orthogonal to the plane (orth) or along the vert.axis of interest (vert.axis).
}
  \item{in.vert.axis}{
This specifies the input vertical projection axis.  
}
  \item{out.vert.axis}{
This specifies the output/requested vertical projection axis.  
}
  
}
\value{
coord.convert: The new coordinates, which will be the same length as the input coord.

beta.convert: The new beta, of length 1.

scat.convert: The new scat, of length 1.

vert.axis.convert: A list with 3 elements. The first is the coord, which will be the same length as the input coord. The second is the new beta, or length 1. The third is the new scat, of length 1.

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

#Here we are assuming our plane forula is z=2x+3y+1 plus Gaussian intrinsic scatter along
#z with sd=4:

excoord.alpha=c(2,3)
exbeta.vert.axis=1
exscat.vert.axis=4

#We want to find the angle theta between the x-axis/y-axis and the z-axis:

print(coord.convert(excoord.alpha, in.coord.type = "alpha", out.coord.type = "theta"))
excoord.theta=coord.convert(excoord.alpha, in.coord.type = "alpha",
out.coord.type = "theta")

#The examples below will give the same answer whether we define our coordinate system in
#alpha or theta, as long as we use the correct 'in.coord.type' option:

print(beta.convert(exbeta.vert.axis, excoord.alpha, in.proj.type='vert.axis',
out.proj.type='orth', in.coord.type='alpha'))
print(beta.convert(exbeta.vert.axis, excoord.theta, in.proj.type='vert.axis',
out.proj.type='orth', in.coord.type='theta'))

#In the above example the distance to the line is negative because the normal vector to
#the plane is pointing away from the origin, so we must move backwards to get there. The
#absolute value of this is the minimum distance between our plane and the origin.

print(scat.convert(exscat.vert.axis, excoord.alpha, in.proj.type='vert.axis',
out.proj.type='orth', in.coord.type='alpha'))
print(scat.convert(exscat.vert.axis, excoord.theta, in.proj.type='vert.axis',
out.proj.type='orth', in.coord.type='theta'))

#In the above example the scatter is less orthogonal to the plane than along the z-axis.
#The scatter orthogonal to the plane will *always* be smaller than the scatter along any
#given axis, so this is what we should expect.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ convert }
\keyword{ hyper }% __ONLY ONE__ keyword per line
