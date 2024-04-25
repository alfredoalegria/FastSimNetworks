dilution_f<-function(x,a){exp(-a^2*x^2)*(2/pi)^0.25*sqrt(a)}

SimDilution2<-function(L,nPointsE,nCopies,scaleParam){
  cholILM=chol(solve(laplacian(L$m,L$dpath)));  
  sites=genSites(nPointsE);
  nsites=nrow(sites);
  zSim=rep(0,nsites);
  for(m in 1:nCopies){
    eps=sample(c(-1,1),1,replace=T);
    x=rcauchy(1,location=0,scale=1);
    zAux=SimAux(lengths_psp(L$lines),nsegments(L),nPointsE,nsites,nrow(cholILM),cholILM,L$from,L$to,sites)[[4]];
    zSim=zSim+eps/sqrt(dcauchy(x))*dilution_f(zAux-x,scaleParam);
  } 
  return(data.frame(e=sites[,2],t=sites[,1],values=zSim*sqrt(1/nCopies)));  
}

