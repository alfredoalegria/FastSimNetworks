library("stlnpp")
library("sp")
dyn.load("/basic_functions/c/BB.so")
source("/basic_functions/r/BB.R")
source("/basic_functions/functions.R")
source("/basic_functions/plotSim.R")
source("/sim_algorithms/auxiliary_process.R")
source("/sim_algorithms/SpectralSim.R")
source("/sim_algorithms/Dilution1Sim.R")
source("/sim_algorithms/Dilution2Sim.R")


set.seed(23);

## Chicago network
data(chicago);
Net=chicago$domain;

## Number of target points (200 per edge)
nPointsE=rep(200,nsegments(Net));  

## Spectral 
DatSim1=SimSpec(Net,nPointsE,nCopies=1000,scaleParam=0.2);
p1=plotSim(DatSim1$values,Net,nPointsE)$plot;

## Poisson dilution 
DatSim2=SimDilution1(Net,nPointsE,nCopies=1000,scaleParam=0.2);
p2=plotSim(DatSim2$values,Net,nPointsE)$plot;

## Dilution of a random germ
DatSim3=SimDilution2(Net,nPointsE,nCopies=1000,scaleParam=0.2);
p3=plotSim(DatSim3$values,Net,nPointsE)$plot;

