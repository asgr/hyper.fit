plot.hyper.fit=function(x,...){
  if(class(x)!='hyper.fit'){stop('Object must be of type hyper.fit')}
  X=x$X
  covarray=x$covarray
  coord=x$parm.vert.axis[1]
  beta=x$parm.vert.axis[2]
  scat=x$parm.vert.axis[3]
  weights=x$weights
  coord.type='alpha'
  proj.type='vert.axis'
  if(x$dims>3){stop('Default plots only exist for 2d/3d data!')}
  if(x$dims==2){hyper.plot2d(X=X, covarray=covarray, fitobj=x, weights=weights, ...)}
  if(x$dims==3){hyper.plot3d(X=X, covarray=covarray, fitobj=x, weights=weights, ...)}
}

hyper.plot2d=function(X,covarray,vars,fitobj,parm.coord,parm.beta,parm.scat,vert.axis,weights,coord.type='alpha',proj.type='orth',errorscale=1,doellipse=TRUE,sigscale=c(0,4),trans=1,dobar=FALSE,position='topright',...){
  if(missing(X)){stop('You must provide X matrix/ data-frame !')}
  checkX=as.numeric(unlist(X))
  if(any(is.null(checkX) | is.na(as.numeric(checkX)) | is.nan(as.numeric(checkX)) | is.infinite(as.numeric(checkX)))){stop('All elements of X must be real numbers, with no NULL, NA, NAN or infinite values.')}
  X=as.matrix(X)
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2] #Number of dimensions
  if(dims<2){stop('The X matrix must have 2 or more columns (i.e. 2 or more dimensions for fitting)')}
  
  if(missing(vert.axis)){vert.axis=dims}
  if(vert.axis>dims){stop(paste('The vertical dimension (',vert.axis,') cannot be higher than the number of dimensions provided (',dims,')',sep=''))}
  
  if(!missing(fitobj)){
    if(class(fitobj)!='hyper.fit'){stop('fitobj class is of the wrong type!. Must be of type \'hyper.fit\'.')}
    coord=fitobj$parm.vert.axis[1]
    beta=fitobj$parm.vert.axis[2]
    scat=fitobj$parm.vert.axis[3]
    coord.type='alpha'
    proj.type='vert.axis'
  }else{
    coord=parm.coord
    beta=parm.beta
    scat=parm.scat
  }
  
  if(missing(covarray)){
    if(!missing(vars)){
      if(any(dim(vars) != dim(X))){stop('vars matrix must be the same size as X matrix!')}
      covarray=array(0,dim = c(dims,dims,N))
      for(i in 1:dims){covarray[i,i,]=vars[,i]}
    }else{
      covarray=array(0,dim = c(dims,dims,N))
    }
  }
  checkcovarray=as.numeric(unlist(covarray))
  if(any(is.null(checkcovarray) | is.na(as.numeric(checkcovarray)) | is.nan(as.numeric(checkcovarray)) | is.infinite(as.numeric(checkcovarray)))){stop('All elements of covarray must be real numbers, with no NULL, NA, NAN or infinite values.')}
  
  if(missing(weights)){weights=rep(1,N)}
  if(length(weights)==1){weights=rep(weights,N)}
  if(length(weights)!=N){stop(paste('weights vector must be the same length as rows in X (',N,')',sep=''))}
  
  #alphas=coord.convert(coord,in.coord.type=coord.type,out.coord.type='alpha')
  #coord.orth=c(alphas,-1)
  #beta.orth=beta.convert(beta=beta,coord=alphas,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  #scat.orth=scat.convert(scat=scat,coord=alphas,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  #beta.vert=beta.convert(beta=beta,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  #scat.vert=scat.convert(scat=scat,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  
  alphas=coord.convert(coord,in.coord.type=coord.type,out.coord.type='alpha')
  beta.vert=beta.convert(beta=beta,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  scat.vert=scat.convert(scat=scat,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  remapaxis=vert.axis.convert(coord=alphas,beta=beta.vert,scat=scat.vert,in.vert.axis=vert.axis,out.vert.axis=dims,in.proj.type='vert.axis',in.coord.type='alpha')
  alphas=remapaxis$coord
  beta.vert=remapaxis$beta
  scat.vert=remapaxis$scat
  coord.orth=remapaxis$vec.orth
  beta.orth=beta.convert(beta=beta.vert,coord=alphas,in.proj.type='vert.axis',out.proj.type='orth',in.coord.type='alpha')
  scat.orth=scat.convert(scat=scat.vert,coord=alphas,in.proj.type='vert.axis',out.proj.type='orth',in.coord.type='alpha')
  
  sx=sqrt(covarray[1,1,])
  sy=sqrt(covarray[2,2,])
  corxy=covarray[1,2,]/(sx*sy)
  
  sigmas=hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='sig',errorscale=errorscale)
  colscale=hsv(magmap(sigmas,lo=sigscale[1],hi=sigscale[2],flip=TRUE,type='num')$map,alpha=trans)
  magplot(X,col=colscale,xlab=colnames(X)[1],ylab=colnames(X)[2],...)
  if(doellipse){magerr(x=X[,1],y=X[,2],xlo=sx,ylo=sy,corxy=corxy,col=colscale,fill=TRUE)}
  abline(beta.vert,alphas)
  abline(beta.vert+scat.vert,alphas,lty=2)
  abline(beta.vert-scat.vert,alphas,lty=2)
  if(dobar){
    colscale=hsv(magmap(seq(sigscale[1],sigscale[2],len=100),lo=sigscale[1],hi=sigscale[2],flip=TRUE,type='num')$map)
    magbar(position=position,range=sigscale,col=colscale,title='Sigma Offset')
  }
  return=sigmas
}

hyper.plot3d=function(X,covarray,vars,fitobj,parm.coord,parm.beta,parm.scat,vert.axis,weights,coord.type='alpha',proj.type='orth',errorscale=1,doellipse=TRUE,sigscale=c(0,4),trans=1,...){
  if(missing(X)){stop('You must provide X matrix/ data-frame !')}
  checkX=as.numeric(unlist(X))
  if(any(is.null(checkX) | is.na(as.numeric(checkX)) | is.nan(as.numeric(checkX)) | is.infinite(as.numeric(checkX)))){stop('All elements of X must be real numbers, with no NULL, NA, NAN or infinite values.')}
  X=as.matrix(X)
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2] #Number of dimensions
  
  if(missing(vert.axis)){vert.axis=dims}
  if(vert.axis>dims){stop(paste('The vertical dimension (',vert.axis,') cannot be higher than the number of dimensions provided (',dims,')',sep=''))}
  
  if(dims<2){stop('The X matrix must have 2 or more columns (i.e. 2 or more dimensions for fitting)')}
  if(!missing(fitobj)){
    if(class(fitobj)!='hyper.fit'){stop('fitobj class is of the wrong type!. Must be of type \'hyper.fit\'.')}
    coord=fitobj$parm.vert.axis[1:2]
    beta=fitobj$parm.vert.axis[3]
    scat=fitobj$parm.vert.axis[4]
    coord.type='alpha'
    proj.type='vert.axis'
  }else{
    coord=parm.coord
    beta=parm.beta
    scat=parm.scat
  }
  
  if(missing(covarray)){
    if(!missing(vars)){
      if(any(dim(vars) != dim(X))){stop('vars matrix must be the same size as X matrix!')}
      covarray=array(0,dim = c(dims,dims,N))
      for(i in 1:dims){covarray[i,i,]=vars[i,]}
    }else{
      covarray=array(0,dim = c(dims,dims,N))
    }
  }
  checkcovarray=as.numeric(unlist(covarray))
  if(any(is.null(checkcovarray) | is.na(as.numeric(checkcovarray)) | is.nan(as.numeric(checkcovarray)) | is.infinite(as.numeric(checkcovarray)))){stop('All elements of covarray must be real numbers, with no NULL, NA, NAN or infinite values.')}
  
  if(missing(weights)){weights=rep(1,N)}
  if(length(weights)==1){weights=rep(weights,N)}
  if(length(weights)!=N){stop(paste('weights vector must be the same length as rows in X (',N,')',sep=''))}
  
  #alphas=coord.convert(coord,in.coord.type=coord.type,out.coord.type='alpha')
  #coord.orth=c(alphas,-1)
  #beta.orth=beta.convert(beta=beta,coord=alphas,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  #scat.orth=scat.convert(scat=scat,coord=alphas,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  #beta.vert=beta.convert(beta=beta,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  #scat.vert=scat.convert(scat=scat,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  
  alphas=coord.convert(coord,in.coord.type=coord.type,out.coord.type='alpha')
  beta.vert=beta.convert(beta=beta,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  scat.vert=scat.convert(scat=scat,coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  remapaxis=vert.axis.convert(coord=alphas,beta=beta.vert,scat=scat.vert,in.vert.axis=vert.axis,out.vert.axis=dims,in.proj.type='vert.axis',in.coord.type='alpha')
  alphas=remapaxis$coord
  beta.vert=remapaxis$beta
  scat.vert=remapaxis$scat
  coord.orth=remapaxis$vec.orth
  beta.orth=beta.convert(beta=beta.vert,coord=alphas,in.proj.type='vert.axis',out.proj.type='orth',in.coord.type='alpha')
  scat.orth=scat.convert(scat=scat.vert,coord=alphas,in.proj.type='vert.axis',out.proj.type='orth',in.coord.type='alpha')
  
  sigmas=hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='sig',errorscale=errorscale)
  colscale=hsv(magmap(sigmas,lo=sigscale[1],hi=sigscale[2],flip=TRUE,type='num')$map,alpha=trans)
  
  plot3d(X,col=colscale,...)
  limx=c(min(X[,1]),max(abs(X[,1])))
  limy=c(min(X[,2]),max(abs(X[,2])))
  limz=c(min(X[,3]),max(abs(X[,3])))
  maxrange=max(c(range(limx),range(limy),range(limz)))
  decorate3d(xlim=c(limx,limx+maxrange),ylim=c(limy,limy+maxrange),zlim=c(limz,limz+maxrange),aspect=1)
  planes3d(alphas[1], alphas[2],-1, beta.vert,alpha=0.2)
  #planes3d(alphas[1], alphas[2],-1, beta.vert+scat.vert,alpha=0.05) #Seems to create glitches, currently removed.
  #planes3d(alphas[1], alphas[2],-1, beta.vert-scat.vert,alpha=0.05) #Seems to create glitches, currently removed.
  
  if(doellipse){
    par3d(skipRedraw=T)
    for(i in 1:length(X[,1])){
      plot3d(ellipse3d(covarray[,,i],centre=X[i,],level=0.6826895,subdivide=2,smooth=T),col=colscale[i],add=T,lit=FALSE,alpha=trans)
    }
    par3d(skipRedraw=F)
  }
  return=sigmas
}
