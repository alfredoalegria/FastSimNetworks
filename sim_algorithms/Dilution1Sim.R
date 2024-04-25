dilution_f<-function(x,a){exp(-a^2*x^2)*(2/pi)^0.25*sqrt(a)}

SimDilution1<-function(L,nPointsE,nCopies,scaleParam){
  cholILM=chol(solve(laplacian(L$m,L$dpath)));  
  sites=genSites(nPointsE);
  nsites=nrow(sites);
  zSim=rep(0,nsites);
  Interv=c(-50,50);
  
  for(m in 1:nCopies){
	  
      zAux=SimAux(lengths_psp(L$lines),nsegments(L),nPointsE,nsites,nrow(cholILM),cholILM,L$from,L$to,sites)[[4]];
      N=rpois(1,lambda=Interv[2]-Interv[1]);
      if(N>0){	  
             eps=cbind(sample(c(-1,1),N,replace=T));
             x=runif(N,min=Interv[1],max=Interv[2]);
             aux=t(sapply(zAux,"-",x));
             zSim=zSim+dilution_f(aux,scaleParam)%*%eps;
      }
  } 
  return(data.frame(e=sites[,2],t=sites[,1],values=zSim*sqrt(1/nCopies)));  
}

