hyper.convert=function(parm,coord,beta=0,scat=0,in.coord.type='alpha',out.coord.type,in.scat.type='vert.axis',out.scat.type,in.vert.axis,out.vert.axis){
  
  #Preamble
  
  if(!missing(parm)){
    dims=length(parm)-1
    if(in.coord.type=='normvec'){
      parmoffset=dims
      coord=parm[1:parmoffset]
      scat=parm[parmoffset+1]
    }else{
      parmoffset=dims-1
      coord=parm[1:parmoffset]
      beta=parm[parmoffset+1]
      scat=parm[parmoffset+2]  
    }
  }else{
    if(missing(coord)){stop('Must provide either parm or coord input!')}
    if(in.coord.type=='normvec'){dims=length(coord);parmoffset=dims}
    if(in.coord.type=='alpha'){dims=length(coord)+1;parmoffset=dims-1}
    if(in.coord.type=='theta'){dims=length(coord)+1;parmoffset=dims-1}
  }
  
  
  if(missing(out.coord.type)){out.coord.type=in.coord.type}
  if(missing(out.scat.type)){out.scat.type=in.scat.type}
  if(missing(in.vert.axis)){in.vert.axis=dims}
  if(missing(out.vert.axis)){out.vert.axis=in.vert.axis}
  
  #In stuff
  
  if(in.coord.type=='normvec'){
    beta.orth=sqrt(sum(coord^2))
    normvec.orth=coord
    unitvec.orth=normvec.orth/beta.orth
  }
  if(in.coord.type=='alpha'){
    normvec.orth=rep(0,dims)
    in.dims=(1:dims)[-in.vert.axis]
    normvec.orth[in.dims]=coord
    normvec.orth[in.vert.axis]=-1
    unitvec.orth=normvec.orth/sqrt(sum(normvec.orth^2))
    beta.orth=beta*unitvec.orth[in.vert.axis]
    unitvec.orth=unitvec.orth*sign(beta.orth)
    beta.orth=abs(beta.orth)
  }
  if(in.coord.type=='theta'){
    normvec.orth=rep(0,dims)
    in.dims=(1:dims)[-in.vert.axis]
    normvec.orth[in.dims]=tan(coord*pi/180)
    normvec.orth[in.vert.axis]=-1
    unitvec.orth=normvec.orth/sqrt(sum(normvec.orth^2))
    beta.orth=beta*unitvec.orth[in.vert.axis]
    unitvec.orth=unitvec.orth*sign(beta.orth)
    beta.orth=abs(beta.orth)
  }
  if(in.scat.type=='orth'){
    scat.orth=abs(scat)
  }
  if(in.scat.type=='vert.axis'){
    scat.orth=abs(scat*unitvec.orth[in.vert.axis])
  }
  
  #unitvec.orth=unitvec.orth*sign(beta.orth)
  #normvec.orth=normvec.orth*sign(beta.orth)
  #beta.orth=beta.orth*sign(beta.orth)
  #print(unitvec.orth)
  #print(beta.orth)
  
  #Out stuff
  
  if(out.coord.type=='normvec'){
    out.dims=1:dims
    out.coord=unitvec.orth*beta.orth
    names(out.coord)=paste('normvec',out.dims,sep='')
  }
  if(out.coord.type=='alpha'){
    out.dims=(1:dims)[-out.vert.axis]
    out.coord=-unitvec.orth/unitvec.orth[out.vert.axis]
    out.coord=out.coord[out.dims]
    names(out.coord)=paste('alpha',out.dims,sep='')
  }
  if(out.coord.type=='theta'){
    out.dims=(1:dims)[-out.vert.axis]
    out.coord=-unitvec.orth/unitvec.orth[out.vert.axis]
    out.coord=atan(out.coord[out.dims])*180/pi
    names(out.coord)=paste('theta',out.dims,sep='')
  }
  if(out.coord.type=='normvec'){
    out.beta=abs(beta.orth)
    names(out.beta)='beta.orth'
  }else{
    out.beta=beta.orth/unitvec.orth[out.vert.axis]
    names(out.beta)='beta.vert'
  }
  if(out.scat.type=='orth'){
    out.scat=scat.orth
    names(out.scat)='scat.orth'
  }
  if(out.scat.type=='vert.axis'){
    out.scat=scat.orth/abs(unitvec.orth[out.vert.axis])
    names(out.scat)='scat.vert'
  }
  
  if(out.coord.type=='normvec'){
    parm=c(out.coord,out.scat)
  }else{
    parm=c(out.coord,out.beta,out.scat)
  }
  
  if(out.coord.type=='normvec' & out.scat.type=='orth'){out.vert.axis=NA}
  
  names(unitvec.orth)=paste('unitvec',1:dims,sep='')
  
  out=list(parm=parm,coord=out.coord,beta=out.beta,scat=out.scat,unitvec=unitvec.orth,beta.orth=beta.orth,scat.orth=scat.orth,coord.type=out.coord.type,scat.type=out.scat.type,vert.axis=out.vert.axis,dims=dims)
  class(out)='hyper.plane.param'
  return=out
}

convert=function(x,coord.type='alpha',scat.type='vert.axis',vert.axis,...){UseMethod("convert")}

convert.hyper.plane.param=function(x,coord.type='alpha',scat.type='vert.axis',vert.axis,...){
  parm=x$parm
  in.coord.type=x$coord.type
  in.scat.type=x$scat.type
  in.vert.axis=x$vert.axis
  if(missing(vert.axis)){vert.axis=x$dims}
  return=hyper.convert(parm=parm,in.coord.type=in.coord.type,out.coord.type=coord.type,in.scat.type=in.scat.type,out.scat.type=scat.type,in.vert.axis=in.vert.axis,out.vert.axis=vert.axis)
}
