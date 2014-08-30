coord.convert=function(coord,in.coord.type='alpha',out.coord.type='alpha'){
  if(in.coord.type==out.coord.type){out=coord}
  if(in.coord.type=='alpha' & out.coord.type=='theta'){out=atan(coord)*180/pi}
  if(in.coord.type=='theta' & out.coord.type=='alpha'){out=tan(coord*pi/180)}
  return=out
}

beta.convert=function(beta,coord,in.proj.type='vert.axis',out.proj.type='vert.axis',in.coord.type='alpha'){
  alphas=coord.convert(coord,in.coord.type=in.coord.type,out.coord.type='alpha')
  if(in.proj.type==out.proj.type){out=beta}
  if(in.proj.type=='orth' & out.proj.type=='vert.axis'){out=-beta*sqrt(sum(alphas^2)+1)}
  if(in.proj.type=='vert.axis' & out.proj.type=='orth'){out=-beta/sqrt(sum(alphas^2)+1)}
  return=out
}

scat.convert=function(scat,coord,in.proj.type='vert.axis',out.proj.type='vert.axis',in.coord.type='alpha'){
  alphas=coord.convert(coord,in.coord.type=in.coord.type,out.coord.type='alpha')
  if(in.proj.type==out.proj.type){out=scat}
  if(in.proj.type=='orth' & out.proj.type=='vert.axis'){out=scat*sqrt(sum(alphas^2)+1)}
  if(in.proj.type=='vert.axis' & out.proj.type=='orth'){out=scat/sqrt(sum(alphas^2)+1)}
  return=out
}