coord.convert=function(coord,in.coord.type='alpha',out.coord.type='alpha'){
  if(in.coord.type==out.coord.type){out=coord}
  if(in.coord.type=='alpha' & out.coord.type=='theta'){
    out=atan(coord)*180/pi
    names(out)=sub('alpha','theta',names(coord))
  }
  if(in.coord.type=='theta' & out.coord.type=='alpha'){
    out=tan(coord*pi/180)
    names(out)=sub('theta','alpha',names(coord))
  }
  return=out
}

beta.convert=function(beta,coord,in.proj.type='vert.axis',out.proj.type='vert.axis',in.coord.type='alpha'){
  alphas=coord.convert(coord,in.coord.type=in.coord.type,out.coord.type='alpha')
  if(in.proj.type==out.proj.type){out=beta}
  betasignmult=-prod(sign(coord))
  if(in.proj.type=='orth' & out.proj.type=='vert.axis'){
    out=-beta*sqrt(sum(alphas^2)+1)
    names(out)='beta.vert'
  }
  if(in.proj.type=='vert.axis' & out.proj.type=='orth'){
    out=-beta/sqrt(sum(alphas^2)+1)
    names(out)='beta.orth'
  }
  return=out
}

scat.convert=function(scat,coord,in.proj.type='vert.axis',out.proj.type='vert.axis',in.coord.type='alpha'){
  alphas=coord.convert(coord,in.coord.type=in.coord.type,out.coord.type='alpha')
  if(in.proj.type==out.proj.type){out=scat}
  if(in.proj.type=='orth' & out.proj.type=='vert.axis'){
    out=scat*sqrt(sum(alphas^2)+1)
    names(out)='scat.vert'
  }
  if(in.proj.type=='vert.axis' & out.proj.type=='orth'){
    out=scat/sqrt(sum(alphas^2)+1)
    names(out)='scat.orth'
  }
  return=out
}

vert.axis.convert=function(coord,beta=0,scat=0,in.vert.axis,out.vert.axis,in.proj.type='vert.axis',in.coord.type='alpha'){
  dims=length(coord)+1
  alphas=coord.convert(coord,in.coord.type=in.coord.type,out.coord.type='alpha')
  scat.orth=scat.convert(scat=scat,coord=alphas,in.proj.type=in.proj.type,out.proj.type='orth',in.coord.type='alpha')
  beta.orth=beta.convert(beta=beta,coord=alphas,in.proj.type=in.proj.type,out.proj.type='orth',in.coord.type='alpha')
  if(in.vert.axis==1){fullcoord=c(-1,alphas[1:(dims-1)])}
  if(in.vert.axis>1 & in.vert.axis<dims){fullcoord=c(alphas[1:(in.vert.axis-1)],-1,alphas[in.vert.axis:(dims-1)])}
  if(in.vert.axis==dims){fullcoord=c(alphas[1:(in.vert.axis-1)],-1)}
  fullcoord= -fullcoord/fullcoord[out.vert.axis]
  coord.final= fullcoord[-out.vert.axis]
  coord.final=coord.convert(coord.final,in.coord.type='alpha',out.coord.type=in.coord.type)
  scat.final=scat.convert(scat=scat.orth,coord=coord.final,in.proj.type='orth',out.proj.type=in.proj.type,in.coord.type=in.coord.type)
  beta.orth=beta.orth*sign(fullcoord[in.vert.axis]*fullcoord[out.vert.axis])
  beta.final=beta.convert(beta=beta.orth,coord=coord.final,in.proj.type='orth',out.proj.type=in.proj.type,in.coord.type=in.coord.type)
  
  dimvec=(1:dims)[-out.vert.axis]
  if(in.coord.type=='alpha'){names(coord.final)=paste('alpha',dimvec,sep='')}
  if(in.coord.type=='theta'){names(coord.final)=paste('theta',dimvec,sep='')}
  
  out=list(coord=coord.final,beta=beta.final,scat=scat.final)
  return=out
}
