hyper.like=function(parm,X,covarray,weights=1,errorscale=1,k.vec=FALSE,output='sum'){
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2]  #Number of dimensions
  coord.orth=parm[1:dims]
  scat.orth=parm[dims+1]
  eTce=projcovarray(covarray,coord.orth)  #The component of covariance along our vector coord.orth
  eTe=sum(coord.orth^2) #Sum of squares of normal vec to renormalise properly
  #orthvariance=scat.orth^2+(eTce/eTe)*(errorscale^2)
  orthvariance=scat.orth^2+(eTce/eTe)*(errorscale^2)
  #if(length(k.vec)> 1){X=X+arrayvecmult(covarray,k.vec)}
  if(length(k.vec)> 1){
    #scat.vec=scat.orth/coord.orth
    scat.vec=scat.orth/(coord.orth/sqrt(eTe))
    X=t(t(X)-k.vec*scat.vec^2)#-arrayvecmult(covarray,k.vec)
  }
  originoffset=(X %*% coord.orth)/sqrt(eTe)-sqrt(eTe)
  #old originoffset=(X %*% coord.orth)/sqrt(eTe)-beta.orth
  if(output=='sum' | output=='val'){
    #Here we fully compute the loglikelihood for all N elements   
    loglike=
      -0.5*(                             #0.5 scaling for all subsequent terms
        log(orthvariance)+                #The normalisation of the Gaussian term
          (originoffset^2)/orthvariance   #The Chi-Sq type term
      )
  }
  if(output=='sig'){
    #Here we fully compute the sigma offset for all N elements   
    sigma=sqrt((originoffset^2)/orthvariance)   #The Chi-Sq type term
  }
  if(output=='sum'){out=sum(weights*loglike)}
  if(output=='val'){out=as.numeric(loglike)}
  if(output=='sig'){out=as.numeric(sigma)}
  return=out
}