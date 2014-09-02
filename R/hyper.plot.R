hyper.plot2d=function(X,covarray,vars,fitobj,parm.coord,parm.beta,parm.scat,coord.type='alpha',proj.type='orth',errorscale=1,doellipse=TRUE,clip=0.05,trans=1,...){
  if(missing(X)){stop('You must provide X matrix!')}
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2] #Number of dimensions
  if(dims<2){stop('The X matrix must have 2 or more columns (i.e. 2 or more dimensions for fitting)')}
  if(!missing(fitobj)){
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
      for(i in 1:dims){covarray[i,i,]=vars[i,]}
    }else{
      covarray=array(0,dim = c(dims,dims,N))
    }
  }
  
  alpha=coord.convert(coord,in.coord.type=coord.type,out.coord.type='alpha')
  coord.orth=c(alpha,-1)
  beta.orth=beta.convert(beta=beta,coord=alpha,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  scat.orth=scat.convert(scat=scat,coord=alpha,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  beta.vert=beta.convert(beta=beta,coord=alpha,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  scat.vert=scat.convert(scat=scat,coord=alpha,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  
  sx=sqrt(covarray[1,1,])
  sy=sqrt(covarray[2,2,])
  corxy=covarray[1,2,]/(sx*sy)
  
  LL=hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,output='sig',errorscale=errorscale)
  colscale=hsv(magmap(LL,lo=clip)$map,alpha=trans)
  magplot(X,col=colscale,...)
  if(doellipse){magerr(x=X[,1],y=X[,2],xlo=sx,ylo=sy,corxy=corxy,col=colscale,fill=TRUE)}
  abline(beta.vert,alpha)
  abline(beta.vert+scat.vert,alpha,lty=2)
  abline(beta.vert-scat.vert,alpha,lty=2)
  
  return=LL
}

hyper.plot3d=function(X,covarray,vars,fitobj,parm.coord,parm.beta,parm.scat,coord.type='alpha',proj.type='orth',errorscale=1,doellipse=TRUE,clip=0.05,trans=1,...){
  if(missing(X)){stop('You must provide X matrix!')}
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2] #Number of dimensions
  if(dims<2){stop('The X matrix must have 2 or more columns (i.e. 2 or more dimensions for fitting)')}
  if(!missing(fitobj)){
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
  
  alpha=coord.convert(coord,in.coord.type=coord.type,out.coord.type='alpha')
  coord.orth=c(alpha,-1)
  beta.orth=beta.convert(beta=beta,coord=alpha,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  scat.orth=scat.convert(scat=scat,coord=alpha,in.proj.type=proj.type,out.proj.type='orth',in.coord.type='alpha')
  beta.vert=beta.convert(beta=beta,coord=alpha,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  scat.vert=scat.convert(scat=scat,coord=alpha,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  
  LL=hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,output='sig',errorscale=errorscale)
  colscale=hsv(magmap(LL,lo=clip)$map,alpha=trans)
  
  plot3d(X,col=colscale,...)
  lim=max(abs(X))
  decorate3d(xlim=c(-lim,lim),ylim=c(-lim,lim),zlim=c(-lim,lim),aspect=1)
  planes3d(alpha[1], alpha[2],-1, beta.vert,alpha=0.2)
  if(doellipse){
    par3d(skipRedraw=T)
    for(i in 1:length(X[,1])){
      plot3d(ellipse3d(covarray[,,i],centre=X[i,],level=0.6826895,subdivide=2,smooth=T),col=colscale[i],add=T,lit=FALSE,alpha=trans)
    }
    par3d(skipRedraw=F)
  }
  return=LL
}
