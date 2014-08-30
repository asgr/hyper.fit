\name{Basic Functions}
\alias{makecovarray2d}
\alias{makecovarray3d}
\alias{makecovmat2d}
\alias{makecovmat3d}
\alias{makerotmat2d}
\alias{makerotmat3d}
\alias{rotdata2d}
\alias{rotdata3d}
\alias{rotcovmat}
\alias{ranrotcovmat2d}
\alias{ranrotcovmat3d}
\alias{makeranrotmat}
\alias{makeformula}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to calculate various basic properties important for line fitting.
}
\description{
rotdata2d: Function to generate a 2xN matrix with rotated columns x and y.

rotdata3d: Function to generate a 3xN matrix with rotated columns x,y and z.

makerotmat2d: Function to generate a 2x2 rotation matrix.

makerotmat3d: Function to generate a 3x3 rotation matrix.

rotcovmat: Function to generate a rotated covariance matrix, either 2x2 or 3x3.

ranrotcovmat2d: Function to generate a randomly rotated 2x2 covariance matrix.

ranrotcovmat3d: Function to generate a randomly rotated 3x3 covariance matrix.

makecovarray2d: Function to generate a 2x2 covariance matrix.

makecovarray3d: Function to generate a 3x3 covariance matrix.

makecovmat2d: Function to generate a 2x2xN covariance array.

makecovmat3d:  Function to generate a 3x3xN covariance array.
}
\usage{
%Utility rotation functions:
rotdata2d(x,y,theta)
rotdata3d(x,y,z,theta,dim='z')
makerotmat2d(theta)
makerotmat3d(theta,dim='z')
rotcovmat(covmat,theta,dim='x')
ranrotcovmat2d(covmat)
ranrotcovmat3d(covmat)
%Utility covariance array/matrix functions:
makecovarray2d(sx,sy,corxy)
makecovarray3d(sx,sy,sz,corxy,corxz,coryz)
makecovmat2d(sx,sy,corxy)
makecovmat3d(sx,sy,sz,corxy,corxz,coryz)
}
\arguments{
  \item{x}{
Vector of x data. Should be the same length as 'y'.
}
  \item{y}{
Vector of y data. Should be the same length as 'x'.
}
  \item{z}{
Vector of z data. Should be the same length as 'x'/'y'.
}
  \item{sx}{
Vector of x errors. Should be the same length as 'sy'.
}
  \item{sy}{
Vector of y errors. Should be the same length as 'sx'.
}
  \item{sz}{
Vector of z errors. Should be the same length as 'sx'/'sy'.
}
  \item{corxy}{
Vector of correlation values between 'sx' and 'sy' values. Should be the same length as 'sx'.
}
  \item{corxz}{
Vector of correlation values between 'sx' and 'sz' values. Should be the same length as 'sx'/'sy'/'sz'/'corxy'/'coryz'.
}
  \item{coryz}{
Vector of correlation values between 'sy' and 'sz' values. Should be the same length as 'sx'/'sy'/'sz'/'corxy'/'corxz'.
}
  \item{covmat}{
A 2x2 or 3x3 covariance matrix that the user wishes to rotate.  
}
  \item{theta}{
Angle in degrees for rotation. x -> y = +ve rotation, x -> z = +ve rotation, y -> z = +ve rotation.
}
  \item{dim}{
In 3D this specifies the axis around which the rotation takes place. If dim='x' rotation is in the yz plane and y -> z = +ve rotation. If dim='y' rotation is in the xz plane and x -> z = +ve rotation. If dim='z' rotation is in the xy plane and x -> y = +ve rotation.
}
}
\value{
rotdata2d: A 2xN matrix with rotated columns x and y.

rotdata3d: A 3xN matrix with rotated columns x,y and z.

makerotmat2d: A 2x2 rotation matrix.

makerotmat3d: A 3x3 rotation matrix.

rotcovmat: A rotated covariance matrix, either 2x2 or 3x3.

ranrotcovmat2d: A randomly rotated 2x2 covariance matrix.

ranrotcovmat3d: A randomly rotated 3x3 covariance matrix.

makecovarray2d: A 2x2 covariance matrix.

makecovarray3d: A 3x3 covariance matrix.

makecovmat2d: A 2x2xN covariance array.

makecovmat3d:  A 3x3xN covariance array.
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
extheta=30    #Experiment with changing this
exdim='z'     #Experiment with chaging to 'x' or 'y'
exvecx=1:10   #Experiment with changin this
exvecy=1:10   #Experiment with changin this
exvecz=1:10   #Experiment with changin this

print(cbind(exvecx,exvecy))
print(rotdata2d(exvecx,exvecy,extheta))
print(rotdata3d(exvecx,exvecy,exvecz,extheta,exdim))
print(makerotmat2d(extheta))
print(makerotmat3d(extheta,dim=exdim))

exsx=1        #Experiment with changin this
exsy=2        #Experiment with changin this
exsz=3        #Experiment with changin this
excorxy=0.8   #Experiment with changin this between -1 and 1
excorxz=-0.3  #Experiment with changin this between -1 and 1
excoryz=0.5   #Experiment with changin this between -1 and 1

print(makecovmat2d(exsx,exsy,excorxy))
print(makecovmat3d(exsx,exsy,exsz,excorxy,excorxz,excoryz))
print(makecovarray2d(exsx*1:4,exsy*1:4,excorxy))
print(makecovarray3d(exsx*1:4,exsy*1:4,exsz*1:4,excorxy,excorxz,excoryz))

excovmat2d=makecovmat2d(exsx,exsy,excorxy)
excovmat3d=makecovmat3d(exsx,exsy,exsz,excorxy,excorxz,excoryz)

print(rotcovmat(excovmat2d,extheta))
print(rotcovmat(excovmat3d,extheta,exdim))
print(ranrotcovmat2d(excovmat2d))
print(ranrotcovmat3d(excovmat3d))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ utility }
\keyword{ rotation }
\keyword{ array }
\keyword{ matrix }
\keyword{ covariance }
