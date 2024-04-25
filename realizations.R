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
nPointsE=rep(200,nsegments(Net));

DatSim1=SimSpec(Net,nPointsE,nCopies=1000,scaleParam=0.2);
p1=plotSim(DatSim1$values,Net,nPointsE)$plot;
DatSim2=SimDilution1(Net,nPointsE,nCopies=1000,scaleParam=0.2);
p2=plotSim(DatSim2$values,Net,nPointsE)$plot;
DatSim3=SimDilution2(Net,nPointsE,nCopies=1000,scaleParam=0.2);
p3=plotSim(DatSim3$values,Net,nPointsE)$plot;

