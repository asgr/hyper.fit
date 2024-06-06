\name{hyper.basic}
\alias{hyper.basic}
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
\alias{projX}
\alias{projcovmat}
\alias{projcovarray}
\alias{arrayvecmult}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to calculate various basic properties important for line fitting
}
\description{
rotdata2d: Function to generate a 2xN matrix with rotated columns x and y.

rotdata3d: Function to generate a 3xN matrix with rotated columns x,y and z.

makerotmat2d: Function to generate a 2x2 rotation matrix.

makerotmat3d: Function to generate a 3x3 rotation matrix.

rotcovmat: Function to generate a rotated covariance matrix, either 2x2 or 3x3.

ranrotcovmat2d: Function to generate a randomly rotated 2x2 covariance matrix.

ranrotcovmat3d: Function to generate a randomly rotated 3x3 covariance matrix.

makecovarray2d: Function to generate a  2x2xN covariance array.

makecovarray3d: Function to generate a  3x3xN covariance array.

makecovmat2d: Function to generate a 2x2 covariance matrix.

makecovmat3d:  Function to generate a 3x3 covariance matrix.

projX:  Projection of position (or possibly velocity) elements along a vector (any number of dimensions, but matrices must conform).

projcovmat: Projection of a covariance matrix along a vector (any number of dimensions, but matrices must conform).

projcovarray: Projection of a covariance array (dim x dim x N) along a vector (any number of dimensions, but matrices must conform).

arrayvecmult: Matrix multiply array elements by a vector. To behave sensibly the second dimension of the AxBxC (i.e. B) array should be the same length as the multiplication vector.

}
\usage{
rotdata2d(x, y, theta)

rotdata3d(x, y, z, theta, dim = 'z')

makerotmat2d(theta)

makerotmat3d(theta, dim = 'z')

rotcovmat(covmat, theta, dim = 'x')

ranrotcovmat2d(covmat)

ranrotcovmat3d(covmat)

makecovarray2d(sx, sy, corxy)

makecovarray3d(sx, sy, sz, corxy, corxz, coryz)

makecovmat2d(sx, sy, corxy)

makecovmat3d(sx, sy, sz, corxy, corxz, coryz)

projX(X, projvec)

projcovmat(covmat, projvec)

projcovarray(covarray, projvec)

arrayvecmult(array, vec)
}
\arguments{
  \item{x}{
Vector of x data. Should be the same length as y.
}
  \item{y}{
Vector of y data. Should be the same length as x.
}
  \item{z}{
Vector of z data. Should be the same length as x/y.
}
  \item{X}{
A position matrix with the N (number of data points) rows by d (number of dimensions) columns.
}
  \item{sx}{
Vector of x errors. Should be the same length as sy.
}
  \item{sy}{
Vector of y errors. Should be the same length as sx.
}
  \item{sz}{
Vector of z errors. Should be the same length as sx/sy.
}
  \item{corxy}{
Vector of correlation values between sx and sy values. Should be the same length as sx.
}
  \item{corxz}{
Vector of correlation values between sx and sz values. Should be the same length as sx/sy/sz/corxy/coryz.
}
  \item{coryz}{
Vector of correlation values between sy and sz values. Should be the same length as sx/sy/sz/corxy/corxz.
}
  \item{covmat}{
A dxd (d=dimensions, i.e. 2x2 for 2d or 3x3 for 3d). The makecovmat2d and makecovmat3d are convenience functions that make populating 2x2 and 3x3 matrices easier for a novice user.  
}
  \item{covarray}{
A dxdxN array containing the full covariance (d=dimensions, N=number of dxd matrices in the array stack). The makecovarray2d and makecovarray3d are convenience functions that make populating 2x2xN and 3x3xN arrays easier for a novice user.
}
  \item{theta}{
Angle in degrees for rotation. x -> y = +ve rotation, x -> z = +ve rotation, y -> z = +ve rotation.
}
  \item{dim}{
In 3D this specifies the axis around which the rotation takes place. If dim='x' rotation is in the yz plane and y -> z = +ve rotation. If dim='y' rotation is in the xz plane and x -> z = +ve rotation. If dim='z' rotation is in the xy plane and x -> y = +ve rotation.
}
  \item{projvec}{
The vector defining the desired projection. This does not need to be of length 1, but the user should be aware that unless it is of length 1 then the solution will not be normalised correctly for the unit vector case (if this is desired the input should be projvec=exprojvec/(sqrt(sum(exprojvec^2))), where exprojvec is the desired direction of projection).  
}
  \item{array}{
A AxBxC array.
}
  \item{vec}{
A vector of length B (same length as the second array dimension of array argument).  
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

projX:  Projection of X along vector (length N).

projcovmat: Projection of covariance matrix along vector (length 1).

projcovarray: Projection of covariance array along vector (length N).

arrayvecmult: Matrix multiplication of array stack slices by vector (size dxN).
}
\references{
Robotham, A.S.G., & Obreschkow, D., PASA, in press
}
\author{
Aaron Robotham and Danail Obreschkow
}
\seealso{
  \code{\link{hyper.basic}}, \code{\link{hyper.convert}}, \code{\link{hyper.fit-data}}, \code{\link{hyper.fit}}, \code{\link{hyper.plot}}, \code{\link{hyper.sigcor}}, \code{\link{hyper.summary}}
}
\examples{
extheta=30    #Experiment with changing this
exdim='z'     #Experiment with chaging to 'x' or 'y'
exvecx=1:10   #Experiment with changin this
exvecy=1:10   #Experiment with changin this
exvecz=1:10   #Experiment with changin this

print(cbind(exvecx, exvecy))
print(rotdata2d(exvecx, exvecy, extheta))
print(rotdata3d(exvecx, exvecy, exvecz, extheta, exdim))
print(makerotmat2d(extheta))
print(makerotmat3d(extheta, dim=exdim))

exsx=1        #Experiment with changing this
exsy=2        #Experiment with changing this
exsz=3        #Experiment with changing this
excorxy=0.8   #Experiment with changing this between -1 and 1
excorxz=-0.3  #Experiment with changing this between -1 and 1
excoryz=0.5   #Experiment with changing this between -1 and 1

print(makecovmat2d(exsx, exsy, excorxy))
print(makecovmat3d(exsx, exsy, exsz, excorxy, excorxz, excoryz))
print(makecovarray2d(exsx*1:4, exsy*1:4, excorxy))
print(makecovarray3d(exsx*1:4, exsy*1:4, exsz*1:4, excorxy, excorxz, excoryz))

excovmat2d=makecovmat2d(exsx, exsy, excorxy)
excovmat3d=makecovmat3d(exsx, exsy, exsz, excorxy, excorxz, excoryz)
excovarray2d=makecovarray2d(exsx*1:4, exsy*1:4, excorxy)
excovarray3d=makecovarray3d(exsx*1:4, exsy*1:4, exsz*1:4, excorxy, excorxz, excoryz)

print(rotcovmat(excovmat2d, extheta))
print(rotcovmat(excovmat3d, extheta, exdim))
print(ranrotcovmat2d(excovmat2d))
print(ranrotcovmat3d(excovmat3d))

exprojvec2d=c(1, 2)
exprojvec2d=exprojvec2d/sqrt(sum(exprojvec2d^2))
exprojvec3d=c(1, 2, 3)
exprojvec3d=exprojvec3d/sqrt(sum(exprojvec3d^2))

print(projX(cbind(exvecx, exvecy), exprojvec2d))
print(projX(cbind(exvecx, exvecy, exvecz), exprojvec3d))
print(projcovmat(excovmat2d, exprojvec2d))
print(projcovmat(excovmat3d, exprojvec3d))
print(projcovarray(excovarray2d, exprojvec2d))
print(projcovarray(excovarray3d, exprojvec3d))

#Notice that the first outputs of the 2d/3d projcovarray example correspond to the outputs
#of the 2d/3d projcovmat examples.

#First for comparison:

print(t(matrix(1:9,3) \%*\% 1:3))
print(t(matrix(10:18,3) \%*\% 1:3))
print(t(matrix(19:27,3) \%*\% 1:3))

#And now an array example of the above operations:

print(arrayvecmult(array(1:27,c(3,3,3)),1:3))

#And an example where all array dimensions are different:

print(matrix(1:6,2) \%*\% 1:3)
print(matrix(7:12,2) \%*\% 1:3)
print(matrix(13:18,2) \%*\% 1:3)
print(matrix(19:24,2) \%*\% 1:3)
print(arrayvecmult(array(1:24,c(2,3,4)),1:3))

#Note that the following is not allowed:

try(arrayvecmult(array(1:24,c(3,2,4)),1:3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ utility }
\keyword{ rotation }
\keyword{ array }
\keyword{ matrix }
\keyword{ covariance }
