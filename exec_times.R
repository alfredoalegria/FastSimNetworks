library("stlnpp")
library("sp")
dyn.load("/Users/alfredo/Desktop/codes_final/basic_functions/c/BB.so")
source("/Users/alfredo/Desktop/codes_final/basic_functions/r/BB.R")
source("/Users/alfredo/Desktop/codes_final/basic_functions/functions.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/auxiliary_process.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/SpectralSim.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/Dilution2Sim.R")
source("/Users/alfredo/Desktop/codes_final/sim_algorithms/Dilution1Sim.R")
source("/Users/alfredo/Desktop/codes_final/basic_functions/plotSim.R")
set.seed(23);

data(chicago);
Net=chicago$domain;

k=2^(5:10);
time=c();
for(j in 1:length(k)){
	nPointsE=rep(k[j],nsegments(Net));
      ini=Sys.time();   
          SimSpec(Net,nPointsE,nCopies=1000,scaleParam=0.2);
      fin=Sys.time();
    time[j]=fin-ini;
    print(time[j]);
}




