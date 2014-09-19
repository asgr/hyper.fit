hyper.sigcor=function(N,parmDF){
  bias2unbias=(sqrt((N-parmDF)/2)*gamma((N-parmDF)/2))/gamma((N-parmDF+1)/2)
  samp2pop=sqrt(N/(N-parmDF))
  return=c(bias2unbias=bias2unbias,samp2pop=samp2pop,sampbias2popunbias=bias2unbias*samp2pop)
}
