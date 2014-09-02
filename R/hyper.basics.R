rotdata2d=function(x,y,theta){
  out=makerotmat2d(theta) %*% rbind(x,y)
  return=cbind(out[1,],out[2,])
}

rotdata3d=function(x,y,z,theta,dim='z'){
  out=makerotmat3d(theta,dim=dim) %*% rbind(x,y,z)
  return=cbind(out[1,],out[2,],out[3,])
}

makerotmat2d=function(theta){
  theta=theta*pi/180
  sintheta=sin(theta)
  costheta=cos(theta)
  return=matrix(c(
    costheta    , -sintheta  ,
    sintheta    , costheta
  ),ncol=2,byrow=T)
}

makerotmat3d=function(theta,dim='z'){
  theta=theta*pi/180
  sintheta=sin(theta)
  costheta=cos(theta)
  
  if(dim=='x' | dim==1){
    out=matrix(c(
      1           , 0           , 0          ,
      0           , costheta    , -sintheta  ,
      0           , sintheta    , costheta
    ),ncol=3,byrow=T)
  }
  
  if(dim=='y' | dim==2){
    out=matrix(c(
      costheta    , 0           , -sintheta  ,
      0           , 1           , 0          ,
      sintheta    , 0           , costheta
    ),ncol=3,byrow=T)
  }
  
  if(dim=='z' | dim==3){
    out=matrix(c(
      costheta    , -sintheta   , 0          ,
      sintheta    , costheta    , 0          ,
      0           , 0           , 1
    ),ncol=3,byrow=T)
  }
  
  return=out
}



rotcovmat=function(covmat,theta,dim='x'){
  if(dim(covmat)[1]==2){rotmat=makerotmat2d(theta)}
  if(dim(covmat)[1]==3){rotmat=makerotmat3d(theta,dim=dim)}
  return=rotmat%*%covmat%*%t(rotmat)
}

ranrotcovmat2d=function(covmat){
  rotmat=makerotmat2d(runif(1,0,360))
  return=rotmat %*% covmat %*% t(rotmat)
}

ranrotcovmat3d=function(covmat){
  
  # draw random 3D rotation matrix from and isotropic distribution
  
  makeranunitvec <- function(){
    # draw random 3D unit vector from an isotropic distribution
    r = runif(1)^(1/3) # radial coordinate
    l = 2*pi*runif(1)  # longitude
    a = 2*runif(1)-1   # cos(latitude)
    return(r*c(cos(l)*sqrt(1-a^2),sin(l)*sqrt(1-a^2),a))
  }
  
  x1 = makeranunitvec()
  x2 = makeranunitvec()
  x3 = 
    c(x1[2]*x2[3]-x1[3]*x2[2],x1[3]*x2[1]-x1[1]*x2[3],x1[1]*x2[2]-x1[2]*x2[1])
  x1 = 
    c(x2[2]*x3[3]-x2[3]*x3[2],x2[3]*x3[1]-x2[1]*x3[3],x2[1]*x3[2]-x2[2]*x3[1])
  rotmat=cbind(x1/sqrt(sum(x1*x1)),x2/sqrt(sum(x2*x2)),x3/sqrt(sum(x3*x3)))
  return=rotmat %*% covmat %*% t(rotmat)
}

makecovarray2d=function(sx,sy,corxy){
  return=aperm(array(c(
    sx^2        , sx*sy*corxy,
    sx*sy*corxy , sy^2
  ),dim = c(length(sx),2,2)),c(2,3,1))
}

makecovarray3d=function(sx,sy,sz,corxy,corxz,coryz){
  return=aperm(array(c(
    sx^2        , sx*sy*corxy , sx*sz*corxz,
    sx*sy*corxy , sy^2        , sy*sz*coryz,
    sx*sz*corxz , sy*sz*coryz , sz^2
  ),dim = c(length(sx),3,3)),c(2,3,1))
}

makecovmat2d=function(sx,sy,corxy){makecovarray2d(sx,sy,corxy)[,,1]}

makecovmat3d=function(sx,sy,sz,corxy,corxz,coryz){makecovarray3d(sx,sy,sz,corxy,corxz,coryz)[,,1]}

projX=function(X,projvec){
  return=abs(as.numeric(rbind(projvec) %*% t(X)))
}

projcovmat=function(covmat,projvec){
  return=abs(as.numeric(rbind(projvec) %*% (covmat %*% cbind(projvec))))
}

projcovarray=function(covarray,projvec){
  #This is the same, but faster than the following using the tensor package:
  #projvec=array(c(1,2),c(length(projvec),1,1))
  #as.numeric(tensor(projvec,tensor(covarray,projvec,2,1),1,1))
  return=abs(as.numeric(colSums(projvec*colSums(aperm(covarray,perm = c(2,1,3))*projvec))))
}

arrayvecmult=function(array,vec){
  if(dim(array)[2] != length(vec)){stop('Non-conformable arguments. Second dimension of array must be the same length as vec.')}
  t(colSums(aperm(array,perm = c(2,1,3))*vec))
}