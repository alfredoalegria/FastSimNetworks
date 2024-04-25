
SimSpec<-function(L,nPointsE,nCopies,scaleParam){
  cholILM=chol(solve(laplacian(L$m,L$dpath)));  
  sites=genSites(nPointsE);
  nsites=nrow(sites);
  zSim=rep(0,nsites);
  for(m in 1:nCopies){
    V=scaleParam;
    U=2*pi*runif(1);
    W=runif(1);
    zAux=SimAux(lengths_psp(L$lines),nsegments(L),nPointsE,nsites,nrow(cholILM),cholILM,L$from,L$to,sites)[[4]];
    zSim=zSim+sqrt(-log(W))*cos(V*zAux+U);
  } 
  return(data.frame(e=sites[,2],t=sites[,1],values=zSim*sqrt(2/nCopies)));  
}


