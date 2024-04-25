genSites <- function(nPointsE){
  e=rep(1:length(nPointsE),nPointsE);
  t=(sequence(nPointsE)-1)/rep(nPointsE-1,nPointsE);
  return(cbind(t,e));
}

laplacian<-function(adjM,vDist){	
   c=adjM/vDist;
   diag(c)=0;
   c0=rowSums(c);
   c0[1]=c0[1]+1;
   return(-c+diag(c0));
}
