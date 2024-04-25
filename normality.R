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

set.seed(2);
data(chicago);
B=owin(c(210,420),c(330,510));
Net=chicago$domain[B];
nPointsE=rep(3,nsegments(Net));
id1=sample(1:sum(nPointsE),2);
id2=sample(1:sum(nPointsE),5);
w=runif(5,min=-10,max=10);
a=0.2;
M=50;

nsamples=100;
nShapiro=100;
level=(1:19)*0.05;

lc1_1=c(); 
lc2_1=c();
lc3_1=c();
lc1_2=c();
lc2_2=c();
lc3_2=c();
res1_1=c();
res2_1=c();
res3_1=c();
res1_2=c();
res2_2=c();
res3_2=c();


for(i in 1:length(level)){
	pp1_1=0;
	pp2_1=0;
	pp3_1=0;
	pp1_2=0;
	pp2_2=0;
	pp3_2=0;
	for(k in 1:nShapiro){
        for(j in 1:nsamples){
             DatSim1=SimSpec(Net,nPointsE,nCopies=M,scaleParam=a);
             DatSim2=SimDilution1(Net,nPointsE,nCopies=M,scaleParam=a);
             DatSim3=SimDilution2(Net,nPointsE,nCopies=M,scaleParam=a);            
             
             lc1_1[j]=sum(w[1:2]*DatSim1[id1,3]);
             lc2_1[j]=sum(w[1:2]*DatSim2[id1,3]);
             lc3_1[j]=sum(w[1:2]*DatSim3[id1,3]); 
             
             lc1_2[j]=sum(w*DatSim1[id2,3]);
             lc2_2[j]=sum(w*DatSim2[id2,3]);
             lc3_2[j]=sum(w*DatSim3[id2,3]);        
        }
        
        pp1_1=pp1_1+ifelse(shapiro.test(lc1_1)$p.value<level[i],1,0);
        pp2_1=pp2_1+ifelse(shapiro.test(lc2_1)$p.value<level[i],1,0);
        pp3_1=pp3_1+ifelse(shapiro.test(lc3_1)$p.value<level[i],1,0);
        
        pp1_2=pp1_2+ifelse(shapiro.test(lc1_2)$p.value<level[i],1,0);
        pp2_2=pp2_2+ifelse(shapiro.test(lc2_2)$p.value<level[i],1,0);
        pp3_2=pp3_2+ifelse(shapiro.test(lc3_2)$p.value<level[i],1,0);
    } 
    res1_1[i]=pp1_1/nShapiro;  
    res2_1[i]=pp2_1/nShapiro;  
    res3_1[i]=pp3_1/nShapiro; 
     
    res1_2[i]=pp1_2/nShapiro;  
    res2_2[i]=pp2_2/nShapiro;  
    res3_2[i]=pp3_2/nShapiro; 
    print(i);
}


#####################################################################################
plot(chicago$domain,main="")
segments(B$xrange[1],B$yrange[1],B$xrange[2],B$yrange[1],col=2,lwd=3,lty=2)
segments(B$xrange[1],B$yrange[2],B$xrange[2],B$yrange[2],col=2,lwd=3,lty=2)
segments(B$xrange[1],B$yrange[1],B$xrange[1],B$yrange[2],col=2,lwd=3,lty=2)
segments(B$xrange[2],B$yrange[1],B$xrange[2],B$yrange[2],col=2,lwd=3,lty=2)
#####################################################################################
plot(Net,main="")
box(col = "red",lty=2,lwd=2)    
ss1=genSites(nPointsE)[id1,]
ss2=genSites(nPointsE)[id2,]

for(j in 1:nrow(ss1)){
w=ss1[j,1];
e=ss1[j,2]
points((1-w)*Net$lines$ends[e,1]+w*Net$lines$ends[e,3],
       (1-w)*Net$lines$ends[e,2]+w*Net$lines$ends[e,4],
       main="",col="red",pch=17,lwd=3,cex=1.5)
}

for(j in 1:nrow(ss2)){
w=ss2[j,1];
e=ss2[j,2]
points((1-w)*Net$lines$ends[e,1]+w*Net$lines$ends[e,3],
       (1-w)*Net$lines$ends[e,2]+w*Net$lines$ends[e,4],
       main="",col="blue",pch=15,lwd=3,cex=1.5)
}
