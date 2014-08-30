hyper.like=function(X,covarray,coord.orth,beta.orth=0,scat.orth=1,errorscale=1,output='sum'){
  N=dim(X)[1] #Number of data points
  dims=dim(X)[2]  #Number of dimensions
  K= -(1/2)*log(2*pi) #scaling factor to make the likelihoods proper
  eTce=abs(colSums(coord.orth*colSums(aperm(covarray,perm = c(2,1,3))*coord.orth)))  #The component of covariance along our vector coord.orth
  eTe=sum(coord.orth^2) #Sum of squares of normal vec to renormalise properly
  orthvariance=scat.orth^2+(eTce/eTe)*(errorscale^2)
  originoffset=(X %*% coord.orth)/sqrt(eTe)-beta.orth
  if(output=='sum' | output=='val'){
    #Here we fully compute the loglikelihood for all N elements   
    loglike=  K -                       #The constant term
      0.5*(                             #0.5 scaling for all subsequent terms
        log(orthvariance)+              #The normalisation of the Gaussian term
          (originoffset^2)/orthvariance   #The Chi-Sq type term
      )
  }
  if(output=='sig'){
    #Here we fully compute the loglikelihood for all N elements   
    loglike=  K -                       #The constant term
      0.5*(                             #0.5 scaling for all subsequent terms
          (originoffset^2)/orthvariance   #The Chi-Sq type term
      )
  }
  if(output=='sum'){out=sum(loglike)}
  if(output=='val'){out=as.numeric(loglike)}
  if(output=='sig'){out=as.numeric(loglike)}
  return=out
}