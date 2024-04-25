BB<-function(nEdges,lenEdges,nPointsE,vecNorm){
   
  storage.mode(lenEdges) <- "double" 
  storage.mode(nPointsE) <- "integer" 
  storage.mode(vecNorm) <- "double" 
  num<-sum(nPointsE);
  o <- .C("BB", nEdges=as.integer(nEdges),lenEdges=lenEdges,nPointsE=nPointsE,
                vecNorm=vecNorm,vec=double(num))
  return(o$vec) 

}
