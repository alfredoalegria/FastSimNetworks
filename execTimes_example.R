library("stlnpp")
library("sp")
dyn.load("/basic_functions/c/BB.so")
source("/basic_functions/r/BB.R")
source("/basic_functions/functions.R")
source("/sim_algorithms/auxiliary_process.R")
source("/sim_algorithms/SpectralSim.R")
source("/sim_algorithms/Dilution1Sim.R")
source("/sim_algorithms/Dilution2Sim.R")
set.seed(23);

## Chicago network
data(chicago);
Net=chicago$domain;

## Number of target points
k=2^(5:10);
time=c();   # it will contain the execution times
for(j in 1:length(k)){
      nPointsE=rep(k[j],nsegments(Net));
      ini=Sys.time();   
      SimSpec(Net,nPointsE,nCopies=1000,scaleParam=0.2);
      fin=Sys.time();
      time[j]=fin-ini;
}




