

SimAux<-function(lenEdges,nEdges,nPointsE,nsites,nVert,cholILM,from,to,sites){
   bridges=BB(nEdges,lenEdges,nPointsE,rnorm(sum(nPointsE))); 
   zVert=t(cholILM)%*%rnorm(nVert);   
   zInterp=(1-sites[,1])*zVert[from[sites[,2]]]+sites[,1]*zVert[to[sites[,2]]];
   zAux=bridges+zInterp; 
   return(list(bridges,zVert,zInterp,zAux));
}  
