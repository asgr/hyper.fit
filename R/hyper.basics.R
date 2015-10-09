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

summary.hyper.fit=function(object,...){
  if(class(object)!='hyper.fit'){stop('Object must be of type hyper.fit')}
  cat(paste('Call used was:\n\n'))
  print(object$call)
  cat(paste('\nData supplied was ',object$N,'rows x ',object$dims,'cols.\n\n',sep=''))
  cat(paste('Requested parameters:\n\n'))
  print(object$parm)
  if(class(object$fit)=='optim' | class(object$fit)=='laplace'){
    if(object$parm.covar[1] != 'singular' & object$parm.covar[1] !='unstable'){
      printerrors=sqrt(diag(object$parm.covar))
    }else{
      printerrors=rep(object$parm.covar[1],length(object$parm))
    }
  }
  if(class(object$fit)=='demonoid'){
    if(object$args$doerrorscale){
      printerrors=object$fit$Summary1[1:(object$dims+2),'SD']
    }else{
      printerrors=object$fit$Summary1[1:(object$dims+1),'SD']
    }
  }
  
  names(printerrors)=paste('err_',names(object$parm),sep='')
  cat(paste('\nErrors for parameters:\n\n'))
  print(printerrors)
  
  if(class(object$fit)=='optim' | class(object$fit)=='laplace'){
  cat(paste('\nThe full parameter covariance matrix:\n\n'))
  print(object$parm.covar)
  }
  
  printLL=object$LL$sum
  names(printLL)='LL'
  cat(paste('\nThe sum of the log-likelihood for the best fit parameters:\n\n'))
  print(printLL)
  

  cat(paste('\nUnbiased population estimator for the intrinsic scatter:\n\n'))
  print(object$sigcor)

  if(class(object$fit)=='demonoid'){
    cat(paste('\nProbability of exactly zero intrinsic scatter:\n\n'))
    print(object$zeroscatprob)
  }
  
  cat(paste('\nStandardised parameters vertically projected along dimension ',object$dims,' (',colnames(object$X)[object$dims],'):\n\n',sep=''))
  print(object$parm.vert.axis)
  parmname=colnames(object$X)
  alphas=object$parm.vert.axis[1:(object$dims-1)]
  beta=object$parm.vert.axis[object$dims]
  scat=object$parm.vert.axis[object$dims+1]
  if(sign(beta)==1){signbeta=' + '}else{signbeta=' - '}
  cat(paste('\nStandardised generative hyperplane equation with unbiased population\nestimator for the intrinsic scatter, vertically projected along dimension ',object$dims,' (',colnames(object$X)[object$dims],'):\n\n',sep=''))
  cat(paste(parmname[object$dims],' ~ N(mu= ',paste(paste(paste(format(alphas,digits=4),parmname[1:(object$dims-1)],sep='*'),collapse=' + ',sep=''),sep=''),signbeta,format(abs(beta),digits=4),' , sigma= ',format(scat*object$sigcor/object$parm[object$dims+1],digits=4),')\n',sep=''))
}
