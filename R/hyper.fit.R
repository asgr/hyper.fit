hyper.fit=function(X,covarray,vars,parm.coord,parm.beta,parm.scat,vert.axis,weights,k.vec,itermax=1e4,coord.type='theta',proj.type='orth',algo.func='optim',algo.method='default',Specs=list(alpha.star=0.44),doerrorscale=FALSE){
  call=match.call(expand.dots = FALSE)
  if(missing(X)){stop('You must provide X matrix/ data-frame !')}
  checkX=as.numeric(unlist(X))
  if(any(is.null(checkX) | is.na(as.numeric(checkX)) | is.nan(as.numeric(checkX)) | is.infinite(as.numeric(checkX)))){stop('All elements of X must be real numbers, with no NULL, NA, NAN or infinite values.')}
  X=as.matrix(X)
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2] #Number of dimensions
  if(dims<2){stop('The X matrix / data-frame must have 2 or more columns (i.e. 2 or more dimensions for fitting)')}
  if(is.null(colnames(X))){colnames(X)=paste('x',1:dims,sep='')}
  
  if(missing(vert.axis)){vert.axis=dims}
  if(vert.axis>dims){stop(paste('The vertical dimension (',vert.axis,') cannot be higher than the number of dimensions provided (',dims,')',sep=''))}
  #Sort out provided starting paramters:
  
  if(missing(covarray)){
    if(!missing(vars)){
      if(any(dim(vars) != dim(X))){stop('vars matrix / data-frame must be the same size as X matrix / data-frame !')}
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
  
  if(missing(parm.coord) & missing(parm.beta) & missing(parm.scat)){
    
    makeformula=function(objname='X',dims,vert.axis,env=.GlobalEnv){
      usedims=(1:dims)[-vert.axis]
      text=paste(objname,'[,',vert.axis,']~',paste(paste(objname,'[,',usedims,']',sep=''),collapse ='+'),sep='')
      return=list(text=text,form=as.formula(text,env=env))
    }
    
    startfit=lm(makeformula('X',dims=dims,vert.axis=vert.axis,env=environment())$form)
    start.beta.vert=startfit$coef[1]
    start.alphas=startfit$coef[2:dims]
    start.scat.vert=sd(startfit$residuals)
    parm.beta=beta.convert(beta=start.beta.vert,coord=start.alphas,in.proj.type='vert.axis',out.proj.type=proj.type,in.coord.type='alpha')
    parm.scat=scat.convert(scat=start.scat.vert,coord=start.alphas,in.proj.type='vert.axis',out.proj.type=proj.type,in.coord.type='alpha')
    parm.coord=coord.convert(start.alphas,in.coord.type='alpha',out.coord.type=coord.type)
  }
  
  if((!missing(parm.coord) & !missing(parm.beta) & !missing(parm.scat))==FALSE){stop('You must provide starting values for all parameters (parm.coord, parm.beta, parm.scat) or none at all.')}
  
  if(length(parm.coord) != dims-1){stop(paste('parm.coord is length',length(parm.coord),'but needs to be length',dims-1))}
  if(length(parm.beta) != 1){stop(paste('parm.beta is length',length(parm.beta),'but needs to be length',1))}
  if(length(parm.scat) != 1){stop(paste('parm.scat is length',length(parm.scat),'but needs to be length',1))}
  if(doerrorscale){
    parm=c(parm.coord,parm.beta,parm.scat,1)
  }else{
    parm=c(parm.coord,parm.beta,parm.scat)
  }
  
  dimvec=(1:dims)[-vert.axis]
  if(coord.type=='alpha' & proj.type=='orth'){
    parm.names=c(paste('alpha',dimvec,sep=''),'beta.orth','scat.orth')
    if(doerrorscale){parm.names=c(parm.names,'errorscale')}
    Data=list(data=list(X=X,covarray=covarray),mon.names='',parm.names=parm.names,N=N,weights=weights,doerrorscale=doerrorscale)
  }
  if(coord.type=='alpha' & proj.type=='vert.axis'){
    parm.names=c(paste('alpha',dimvec,sep=''),'beta.vert','scat.vert')
    if(doerrorscale){parm.names=c(parm.names,'errorscale')}
    Data=list(data=list(X=X,covarray=covarray),mon.names='',parm.names=parm.names,N=N,weights=weights,doerrorscale=doerrorscale)
  }
  if(coord.type=='theta' & proj.type=='orth'){
    parm.names=c(paste('theta',dimvec,sep=''),'beta.orth','scat.orth')
    if(doerrorscale){parm.names=c(parm.names,'errorscale')}
    Data=list(data=list(X=X,covarray=covarray),mon.names='',parm.names=parm.names,N=N,weights=weights,doerrorscale=doerrorscale)
  }
  if(coord.type=='theta' & proj.type=='vert.axis'){
    parm.names=c(paste('theta',dimvec,sep=''),'beta.vert','scat.vert')
    if(doerrorscale){parm.names=c(parm.names,'errorscale')}
    Data=list(data=list(X=X,covarray=covarray),mon.names='',parm.names=parm.names,N=N,weights=weights,doerrorscale=doerrorscale)
  }

  if(!missing(k.vec)){
    if(length(k.vec) != dims){stop(paste('The length of k.vec (',length(k.vec),') must be the same as the dimensions provided (',dims,')',sep=''))}
    Data$k.vec=k.vec
  }else{Data$k.vec=FALSE;k.vec=NULL}
  
  argoptions=list(vert.axis=vert.axis,weights=weights,k.vec=k.vec,itermax=itermax,coord.type=coord.type,proj.type=proj.type,algo.func=algo.func,algo.method=algo.method,Specs=Specs,doerrorscale=doerrorscale)
  
  linelikemodel=function(parm,Data){
    if(vert.axis==1){coord.orth=c(-1,parm[1:(dims-1)])}
    if(vert.axis>1 & vert.axis<dims){coord.orth=c(parm[1:(vert.axis-1)],-1,parm[vert.axis:(dims-1)])}
    if(vert.axis==dims){coord.orth=c(parm[1:(vert.axis-1)],-1)}
    if(coord.type=='theta'){coord.orth=tan(coord.orth*pi/180);coord.orth[vert.axis]=-1}
    beta=parm[dims]
    if(proj.type=='vert.axis'){beta.orth= -beta/sqrt(sum(coord.orth^2))}else{beta.orth=beta}
    scat=parm[dims+1]
    if(scat<0){scat=0;parm[dims+1]=0}
    if(Data$doerrorscale){
      parm[dims+2]=abs(parm[dims+2])
      errorscale=parm[dims+2]
    }else{errorscale=1}
    if(proj.type=='vert.axis'){scat.orth=scat/sqrt(sum(coord.orth^2))}else{scat.orth=scat}
    LL=hyper.like(X=Data$data$X, covarray=Data$data$covarray, coord.orth=coord.orth, beta.orth=beta.orth, scat.orth=scat.orth, weights=Data$weights, errorscale=errorscale, k.vec=Data$k.vec, output='sum')
    LP=LL
    if(algo.func=='optim'){out=LP}
    if(algo.func=='LA' | algo.func=='LD'){out=list(LP=LP,Dev=2*LL,Monitor=1,yhat=1,parm=parm)}
    return=out
  }
  #Run the optimisers
  if(algo.func=='optim'){
    if(algo.method=='default'){algo.method='Nelder-Mead'}
    fit=optim(par=parm,fn=linelikemodel,Data=Data,hessian=TRUE,control=list(fnscale=-1),method=algo.method)
    hessdims=dim(fit$hessian)[1]
    fit$hessian[rbind((1:hessdims)[-(dims+1)],(1:hessdims)[-(dims+1)])]=fit$hessian[rbind((1:hessdims)[-(dims+1)],(1:hessdims)[-(dims+1)])]*sign(parm[dims+1])
    fit$Covar=try(solve(fit$hessian))
    if(class(fit$Covar) != 'try-error'){fit$Covar=fit$Covar*sign(fit$Covar[1,1])}
    parm=fit$par
    names(parm)=parm.names
    parm[dims+1]=abs(parm[dims+1])
  }
  if(algo.func=='LA'){
    if(algo.method=='default'){algo.method='NM'}
    fit=LaplaceApproximation(linelikemodel,Data=Data,Iterations=itermax,Method=algo.method,parm=parm,CovEst="Hessian")
    if(doerrorscale){getelements=dims+2}else{getelements=dims+1}
    parm=fit$Summary1[1:getelements,'Mode']
  }
  if(algo.func=='LD'){
    if(algo.method=='default'){algo.method='CHARM'}
    fit=LaplacesDemon(linelikemodel,Data=Data,Iterations=itermax,Algorithm=algo.method,Initial.Values=parm,Specs=Specs,Status=itermax/10)
    if(doerrorscale){getelements=dims+2}else{getelements=dims+1}
    parm=fit$Summary1[1:getelements,'Mean']
    zeroscatprob=length(which(fit$Posterior1[,dims+1]<=0))/length(fit$Posterior1[,1])
  }
  #Prepare output data
  alphas=coord.convert(parm[1:(dims-1)],in.coord.type=coord.type,out.coord.type='alpha')
  beta.vert=beta.convert(beta=parm[dims],coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  scat.vert=scat.convert(scat=parm[dims+1],coord=alphas,in.proj.type=proj.type,out.proj.type='vert.axis',in.coord.type='alpha')
  remapaxis=vert.axis.convert(coord=alphas,beta=beta.vert,scat=scat.vert,in.vert.axis=vert.axis,out.vert.axis=dims,in.proj.type='vert.axis',in.coord.type='alpha')
  alphas=remapaxis$coord
  beta.vert=remapaxis$beta
  scat.vert=remapaxis$scat
  coord.orth=c(alphas,-1)
  beta.orth=beta.convert(beta=beta.vert,coord=alphas,in.proj.type='vert.axis',out.proj.type='orth',in.coord.type='alpha')
  scat.orth=scat.convert(scat=scat.vert,coord=alphas,in.proj.type='vert.axis',out.proj.type='orth',in.coord.type='alpha')
  
  if(doerrorscale){
    LLsum= hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='sum',errorscale=parm[dims+2])
    LLval= hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='val',errorscale=parm[dims+2])
    LLsig= -2*hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='sig',errorscale=parm[dims+2])
  }else{
    LLsum= hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='sum',errorscale=1)
    LLval= hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='val',errorscale=1)
    LLsig= hyper.like(X=X,covarray=covarray,coord.orth=coord.orth,beta.orth=beta.orth,scat.orth=scat.orth,weights=weights,output='sig',errorscale=1)
  }
  
  
  if(doerrorscale){
    parm.vert.axis=c(alphas,beta.vert,scat.vert,parm[dims+2])
    names(parm.vert.axis)=c(paste('alpha',1:(dims-1),sep=''),'beta.vert','scat.vert','errorscale')
  }else{
    parm.vert.axis=c(alphas,beta.vert,scat.vert)
    names(parm.vert.axis)=c(paste('alpha',1:(dims-1),sep=''),'beta.vert','scat.vert')
  }
  
    
  if(algo.func=='optim' | algo.func=='LA'){
    if(class(fit$Covar) != 'try-error'){
      parm.covar=fit$Covar
      colnames(parm.covar)=names(parm)
      rownames(parm.covar)=names(parm)
    }else{parm.covar='singular'}
    sigcor=hyper.sigcor(N,length(parm)-1)[3]*as.numeric(parm[dims+1])
    names(sigcor)=names(parm[dims+1])
    out=list(parm=parm,parm.vert.axis=parm.vert.axis,parm.covar=parm.covar,fit=fit,sigcor=sigcor)
  }
  if(algo.func=='LD'){
    out=list(parm=parm,parm.vert.axis=parm.vert.axis,fit=fit,zeroscatprob=zeroscatprob)
  }
  if(doerrorscale){
    out$parm['errorscale']=abs(out$parm['errorscale'])
    out$parm.vert.axis['errorscale']=abs(out$parm.vert.axis['errorscale'])
  }
  out$N=N
  out$dims=dims
  out$X=X
  out$covarray=covarray
  out$weights=weights
  out$call=call
  out$args=argoptions  
  out$LL=list(sum=LLsum,val=LLval,sig=LLsig)
  
  class(out)='hyper.fit'
  if(algo.func=='optim'){class(out$fit)='optim'} 
  return=out
}
