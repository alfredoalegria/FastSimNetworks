library("stlnpp")
library("sp")
dyn.load("/Users/alfredo/Desktop/codes_final/basic_functions/c/BB.so")
source("/Users/alfredo/Desktop/codes_final/basic_functions/r/BB.R")
dyn.load("/Users/alfredo/Desktop/codes_final/basic_functions/c/resistance_metric.so")
source("/Users/alfredo/Desktop/codes_final/basic_functions/r/resistance_metric.R")
source("/Users/alfredo/Desktop/codes_final/basic_functions/functions.R")
source("/Users/alfredo/Desktop/codes_final/basic_functions/empVariog.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/auxiliary_process.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/SpectralSim.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/Dilution2Sim.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/Dilution1Sim.R")
source("/Users/alfredo/Desktop/codes_final/basic_functions/plotSim.R")
set.seed(23);

## Chicago network
data(chicago)
Net=chicago$domain;

## Generate locations (4 per edge) and select interior points (2 per edge)
nPointsE=rep(4,nsegments(Net));
sites=genSites(nPointsE);
q=sites[,1]*(1-sites[,1])!=0;
sites=sites[q,];

## Laplacian matrix
invL=solve(laplacian(Net$m,Net$dpath));

## Matrix of distances
distMat=resistance_metric(nvertices(Net),nrow(sites),sites[,2],sites[,1],
                          c(invL),Net$from,Net$to,lengths_psp(Net$lines));

## Parameters
scaleParam=0.2;  # correlation parameter
nCopies=1000;    # M value in the central limit approximation
nlags=26;        # number of lags for empirical variogram
hmax=250;        # maximum distance for empirical variogram
tol=0.8;         # distance tolerance for empirical variogram
nrep=200;        # number of realizations
algorithm=1;     # 1=spectral; 2=poisson_dilution; 3=germ_dilution

lags=seq(0,hmax,l=nlags);
vv=array(dim=c(nrep,nlags));  # it contains the empirical variograms
mm=array(dim=c(nrep,nlags));  # it contains the empirical madograms

for(rep in 1:nrep){
   if(algorithm==1){zVec=SimSpec(Net,nPointsE,nCopies,scaleParam)$values;}
   if(algorithm==2){zVec=SimDilution1(Net,nPointsE,nCopies,scaleParam)$values;}
   if(algorithm==3){zVec=SimDilution2(Net,nPointsE,nCopies,scaleParam)$values;}
   vv[rep,]=empVar(nlags,hmax,tol,zVec[q],distMat)[1,];
   mm[rep,]=empVar(nlags,hmax,tol,zVec[q],distMat)[2,];
}

## T values of the statistical test 
t_value=c();  # t values for variograms
t_value2=c(); # t values for madograms
id=c(2,6,11,16,21,26); # lags considered (10,50,100,150,200,250)

for(j in 1:length(id)){
	if(algorithm==1){aux=1-exp(-scaleParam^2*lags[id[j]]/2)}
	if(algorithm>1){aux=1-1/sqrt(1+scaleParam^2*lags[id[j]])}
        t_value[j]= sqrt(199)*(mean(0.5*vv[,id[j]])-aux)/sd(0.5*vv[,id[j]]);
	t_value2[j]= sqrt(199)*(mean(0.5*mm[,id[j]])-sqrt(aux/pi))/sd(0.5*mm[,id[j]]);
}
