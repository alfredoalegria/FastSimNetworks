resistance_metric<-function(nVert,nsites,sites_e,sites_t,ILM,from,to,lenEdges){
   
  storage.mode(sites_e) <- "integer"
  storage.mode(sites_t) <- "double"
  storage.mode(ILM) <- "double"
  storage.mode(from) <- "integer"
  storage.mode(to) <- "integer"
  storage.mode(lenEdges) <- "double" 

  num<-nsites^2;
  o <- .C("resistance_metric", nVert=as.integer(nVert),nsites=as.integer(nsites),sites_e=sites_e,sites_t=sites_t,
                               ILM=ILM,from=from,to=to,lenEdges=lenEdges,vec=double(num))
  return(array(dim=c(nsites,nsites),o$vec)) 

}
