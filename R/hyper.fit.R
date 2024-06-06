hyper.fit=function(X,covarray,vars,parm,parm.coord,parm.beta,parm.scat,parm.errorscale=1,vert.axis,weights,k.vec,prior,itermax=1e4,coord.type='alpha',scat.type='vert.axis',algo.func='optim',algo.method='default',Specs=list(Grid=seq(-0.1,0.1, len=5),dparm=NULL, CPUs=1, Packages=NULL, Dyn.libs=NULL),doerrorscale=FALSE,...){
  #MEGA TEDIOUS CODE PREAMBLE (CATCHING ERRORS AND INITIALISING STUFF)
  
  call=match.call(expand.dots = FALSE)
  if(missing(X)){stop('You must provide X matrix / data-frame !')}
  checkX=as.numeric(unlist(X))
  if(any(is.null(checkX) | is.na(as.numeric(checkX)) | is.nan(as.numeric(checkX)) | is.infinite(as.numeric(checkX)))){stop('All elements of X must be real numbers, with no NULL, NA, NAN or infinite values.')}
  X=as.matrix(X)
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2] #Number of dimensions
  if(dims<2){stop('The X matrix / data-frame must have 2 or more columns (i.e. 2 or more dimensions for fitting)')}
  if(is.null(colnames(X))){colnames(X)=paste('x',1:dims,sep='')}
  
  if(coord.type=='normvec'){parmoffset=dims-1}
  if(coord.type=='alpha'){parmoffset=dims-1}
  if(coord.type=='theta'){parmoffset=dims-1}
  
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
  
  if(!sum(missing(parm.coord),missing(parm.beta),missing(parm.scat)) %in% c(0,3)){stop('You must provide starting values for all parameters via parm.coord / parm.beta / parm.scat or none at all.')}
  if(!missing(parm) & any(!missing(parm.coord) & !missing(parm.beta) & !missing(parm.scat))){stop('You must not supply parm.coord / parm.beta / parm.scat if using parm')}
  
  if(doerrorscale){parmextra=3}else{parmextra=2}
  
  if(!missing(parm)){
    if(length(parm)!=parmoffset+parmextra){stop(paste('If supplying parm it must be of length dimensions+',parmextra-1,' (',parmoffset+parmextra,'). Currently it is of length ',length(parm),'.',sep=''))}
    start.fit=hyper.convert(parm[1:(parmoffset+2)],in.coord.type=coord.type,out.coord.type=coord.type,in.scat.type=scat.type,out.scat.type=scat.type,in.vert.axis=vert.axis,out.vert.axis=vert.axis)
    if(doerrorscale){parm.errorscale=parm[parmoffset+3]}
  }
  
  if(missing(parm) & missing(parm.coord) & missing(parm.beta) & missing(parm.scat)){
    makeformula=function(objname='X',dims,vert.axis){
      usedims=(1:dims)[-vert.axis]
      text=paste(objname,'[,',vert.axis,']~',paste(paste(objname,'[,',usedims,']',sep=''),collapse ='+'),sep='')
      return=list(text=text,form=as.formula(text))
    }
    if(any(covarray[vert.axis,vert.axis,]>0)){
      startweights=1/sqrt(covarray[vert.axis,vert.axis,])
      startweights[is.infinite(startweights)]=max(startweights[is.finite(startweights)])
    }else{
      startweights=rep(1,N)
    }
    startfit=lm(makeformula('X',dims=dims,vert.axis=vert.axis)$form,weights=startweights*weights)
    startfit$coef[is.na(startfit$coef)]=0
    start.alphas=startfit$coef[2:dims]
    start.beta.vert=startfit$coef[1]
    start.scat.vert=sd(startfit$residuals)
    start.fit=hyper.convert(coord=start.alphas,beta=start.beta.vert,scat=start.scat.vert,in.coord.type='alpha',out.coord.type=coord.type,in.scat.type='vert.axis',out.scat.type=scat.type,in.vert.axis=dims,out.vert.axis=vert.axis)
  }
  
  if(missing(parm) & !missing(parm.coord) & !missing(parm.beta) & !missing(parm.scat)){
    start.fit=hyper.convert(coord=parm.coord,beta=parm.beta,scat=parm.scat,in.coord.type=coord.type,out.coord.type=coord.type,in.scat.type=scat.type,out.scat.type=scat.type,in.vert.axis=vert.axis,out.vert.axis=vert.axis)
  }

  if(doerrorscale){
    names(parm.errorscale)='errorscale'
    parm=c(start.fit$parm,parm.errorscale)
  }else{
    parm=c(start.fit$parm)
  }
  parm.names=names(parm)
  if(missing(prior)){prior=function(x){0}}
  Data=list(data=list(X=X,covarray=covarray),mon.names='',parm.names=parm.names,N=N,weights=weights,doerrorscale=doerrorscale,algo.func=algo.func,prior=prior)

  if(!missing(k.vec)){
    if(length(k.vec) != dims){stop(paste('The length of k.vec (',length(k.vec),') must be the same as the dimensions provided (',dims,')',sep=''))}
    Data$k.vec=k.vec
  }else{Data$k.vec=FALSE;k.vec=NULL}
  
  argoptions=list(vert.axis=vert.axis,weights=weights,k.vec=k.vec,itermax=itermax,coord.type=coord.type,scat.type=scat.type,algo.func=algo.func,algo.method=algo.method,Specs=Specs,doerrorscale=doerrorscale)
  
  # END OF MEGA TEDIOUS CODE PREAMBLE (CATCHING ERRORS AND INITIALISING STUFF)

  linelikemodel=function(parm,Data){
    if(Data$algo.func=='LD'){fixparm=TRUE}else{fixparm=FALSE}
    if(fixparm){
      if(parm[parmoffset+2]<0){parm[parmoffset+2]=0}
    }
    if(Data$doerrorscale){
      errorscale=abs(parm[parmoffset+3])
      if(fixparm){parm[parmoffset+3]=errorscale}
    }else{errorscale=1}
    convert.out=hyper.convert(parm=parm[1:(parmoffset+2)],in.coord.type=coord.type,out.coord.type='normvec',in.scat.type=scat.type,out.scat.type='orth',in.vert.axis=vert.axis,out.vert.axis=vert.axis)$parm
    LL=hyper.like(parm=convert.out, X=Data$data$X, covarray=Data$data$covarray, weights=Data$weights, errorscale=errorscale, k.vec=Data$k.vec, output='sum')
    LP=LL+Data$prior(parm)
    if(algo.func=='optim'){out=LP}
    if(algo.func=='LA' | algo.func=='LD'){out=list(LP=LP,Dev=-2*LL,Monitor=1,yhat=1,parm=parm)}
    return=out
  }
  
  #Run the optimisers
  if(algo.func=='optim'){
    if(algo.method=='default'){algo.method='Nelder-Mead'}
    fit=optim(par=parm,fn=linelikemodel,Data=Data,hessian=TRUE,control=list(fnscale=-1),method=algo.method,...)
    parm=fit$par
    names(parm)=parm.names
    fit$Covar=suppressWarnings(try(solve(fit$hessian)))
    if(!inherits(fit$Covar,'try-error')){
      fit$Covar=fit$Covar*sign(fit$Covar[1,1])
    }
    if(parm[parmoffset+2]<0){
      parm[parmoffset+2]=abs(parm[parmoffset+2])
      fit$par[parmoffset+2]=abs(fit$par[parmoffset+2])
    }
  }

  if(algo.func=='LA'){
    if(algo.method=='default'){algo.method='NM'}
    fit=LaplaceApproximation(linelikemodel,Data=Data,Iterations=itermax,Method=algo.method,parm=parm,CovEst="Hessian",...)
    if(doerrorscale){getelements=parmoffset+3}else{getelements=parmoffset+2}
    parm=fit$Summary1[1:getelements,'Mode']
    if(parm[parmoffset+2]<0){
      parm[parmoffset+2]=abs(parm[parmoffset+2])
      fit$Summary1[parmoffset+2,]=suppressWarnings(try(abs(fit$Summary1[parmoffset+2,])))
      fit$Summary2[parmoffset+2,]=suppressWarnings(try(abs(fit$Summary2[parmoffset+2,])))
      covdims=dim(fit$Covar)[1]
      fit$Covar[parmoffset+2,(1:covdims)[-(parmoffset+2)]]=-fit$Covar[parmoffset+2,(1:covdims)[-(parmoffset+2)]]
      fit$Covar[(1:covdims)[-(parmoffset+2)],parmoffset+2]=-fit$Covar[(1:covdims)[-(parmoffset+2)],parmoffset+2]
    }
  }

  if(algo.func=='LD'){
    if(algo.method=='default'){algo.method='GG'}
    fit=LaplacesDemon(linelikemodel,Data=Data,Iterations=itermax,Algorithm=algo.method,Initial.Values=parm,Specs=Specs,...)
    if(doerrorscale){getelements=parmoffset+3}else{getelements=parmoffset+2}
    parm=fit$Summary1[1:getelements,'Mean']
    if(parm[parmoffset+2]<0){
      parm[parmoffset+2]=abs(parm[parmoffset+2])
      fit$Posterior1[,parmoffset+2]=suppressWarnings(try(abs(fit$Posterior1[,parmoffset+2])))
      fit$Posterior2[,parmoffset+2]=suppressWarnings(try(abs(fit$Posterior2[,parmoffset+2])))
      fit$Summary1[parmoffset+2,]=suppressWarnings(try(abs(fit$Summary1[parmoffset+2,])))
      fit$Summary2[parmoffset+2,]=suppressWarnings(try(abs(fit$Summary2[parmoffset+2,])))
    }
    zeroscatprob=length(which(fit$Posterior1[,parmoffset+2]<=0))/length(fit$Posterior1[,1])
    names(zeroscatprob)='zeroscatprob'
  }

  #Prepare output data
  final.coord=parm[1:parmoffset]
  final.beta=parm[parmoffset+1]
  final.scat=parm[parmoffset+2]
  final.fit=hyper.convert(parm=parm[1:(parmoffset+2)],in.coord.type=coord.type,out.coord.type=coord.type,in.scat.type=scat.type,out.scat.type=scat.type,in.vert.axis=vert.axis,out.vert.axis=vert.axis)
  #final.fit.vert=hyper.convert(parm=parm[1:(parmoffset+2)],in.coord.type=coord.type,out.coord.type='alpha',in.scat.type=scat.type,out.scat.type='vert.axis',in.vert.axis=vert.axis,out.vert.axis=dims)
  final.fit.vert=convert(final.fit,coord.type='alpha',scat.type='vert.axis')
  alphas=final.fit.vert$coord
  beta.vert=final.fit.vert$beta
  scat.vert=final.fit.vert$scat
  #final.fit.normvec.orth=hyper.convert(parm=parm[1:(parmoffset+2)],in.coord.type=coord.type,out.coord.type='normvec',in.scat.type=scat.type,out.scat.type='orth',in.vert.axis=vert.axis,out.vert.axis=vert.axis)
  final.fit.normvec.orth=convert(final.fit,coord.type='normvec',scat.type='orth')
  
  if(doerrorscale){
    LLsum= hyper.like(parm=final.fit.normvec.orth$parm,X=X,covarray=covarray,weights=weights,output='sum',errorscale=parm[parmoffset+3])
    LLval= hyper.like(parm=final.fit.normvec.orth$parm,X=X,covarray=covarray,weights=weights,output='val',errorscale=parm[parmoffset+3])
    LLsig= -2*hyper.like(parm=final.fit.normvec.orth$parm,X=X,covarray=covarray,weights=weights,output='sig',errorscale=parm[parmoffset+3])
  }else{
    LLsum= hyper.like(parm=final.fit.normvec.orth$parm,X=X,covarray=covarray,weights=weights,output='sum',errorscale=1)
    LLval= hyper.like(parm=final.fit.normvec.orth$parm,X=X,covarray=covarray,weights=weights,output='val',errorscale=1)
    LLsig= hyper.like(parm=final.fit.normvec.orth$parm,X=X,covarray=covarray,weights=weights,output='sig',errorscale=1)
  }
  
  
  if(doerrorscale){
    parm.vert.axis=c(alphas,beta.vert,scat.vert,parm[parmoffset+3])
    #names(parm.vert.axis)=c(paste('alpha',1:(dims-1),sep=''),'beta.vert','scat.vert','errorscale')
  }else{
    parm.vert.axis=c(alphas,beta.vert,scat.vert)
    #names(parm.vert.axis)=c(paste('alpha',1:(dims-1),sep=''),'beta.vert','scat.vert')
  }
  
  if(algo.func=='optim' | algo.func=='LA'){
    if(!inherits(fit$Covar,'try-error')){
      parm.covar=fit$Covar
      colnames(parm.covar)=names(parm)
      rownames(parm.covar)=names(parm)
    }else{parm.covar='singular'}
    sigcor=hyper.sigcor(N,length(parm)-1)[3]*as.numeric(parm[parmoffset+2])
    names(sigcor)=names(parm[parmoffset+2])
    out=list(parm=parm,parm.vert.axis=parm.vert.axis,fit=fit,sigcor=sigcor,parm.covar=parm.covar)
  }
  if(algo.func=='LD'){
    sigcor=hyper.sigcor(N,length(parm)-1)[1]*as.numeric(parm[parmoffset+2])
    names(sigcor)=names(parm[parmoffset+2])
    out=list(parm=parm,parm.vert.axis=parm.vert.axis,fit=fit,sigcor=sigcor,zeroscatprob=zeroscatprob)
  }
  if(doerrorscale){
    out$parm['errorscale']=abs(out$parm['errorscale'])
    out$parm.vert.axis['errorscale']=abs(out$parm.vert.axis['errorscale'])
  }
  out$hyper.plane=final.fit
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
  
  out$func=function(x, dim='all', parm.vert.axis){
    output=rep(0, dim(x)[1])
    if(doerrorscale){
      parm.limit = 3
    }else{
      parm.limit = 2
    }
    for(i in 1:(length(parm.vert.axis)-parm.limit)){
      if(dim[1]=='all' | i %in% dim){
        output = output + x[,i]*parm.vert.axis[i]
      }
    }
    output = output + parm.vert.axis[i+1]
    return(output)
  }
  
  formals(out$func)$parm.vert.axis = out$parm.vert.axis
  
  out$predict.vert.axis = out$func(out$X)
  
  return=out
}
