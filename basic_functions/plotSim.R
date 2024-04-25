plotSim<-function(simValues,Net,nPointsE){
  e=Net$lines$ends;
  e.x=c();
  e.y=c();   
  for(j in 1:nrow(e)){
	 w=(sequence(nPointsE[j])-1)/rep(nPointsE[j]-1,nPointsE[j]);
     e.x=c(e.x,(1-w)*e[j,1]+w*e[j,3]);
     e.y=c(e.y,(1-w)*e[j,2]+w*e[j,4]);
  }
  df=data.frame(xx=e.x,yy=e.y,zz=simValues);
  coordinates(df)=~xx+yy;
  spp=spplot(df,col="transparent",colorkey=TRUE,cex=0.1,pch=20);
  return(list(plot=spp,coord=cbind(e.x,e.y)));
}

