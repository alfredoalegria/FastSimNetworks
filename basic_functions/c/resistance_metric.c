#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

void resistance_metric(int *nVert,int *nsites,int *sites_e,double *sites_t,double *ILM,int *from,int *to,double *lenEdges,double *vec)
{  
     int id=0,segu,segv;
     double aux,tu,tv,li,lj;
     double Sigajaj,Sigajbj,Sigbjbj,Sigaiai,Sigaibi,Sigbibi,Sigajai,Sigajbi,Sigbjai,Sigbjbi;
  
     for(int i=0;i<(*nsites);i++)
     {  
         for(int j=0;j<(*nsites);j++)
         {    
			  if(i==j){vec[id]=0;}
			  else{			         
			           
			           segu = sites_e[i]-1;
			           segv = sites_e[j]-1;
			           tu = sites_t[i];
			           tv = sites_t[j];
			           if(segu==segv){aux=1;}
			           if(segu!=segv){aux=0;}
			         
			           Sigajaj = ILM[(from[segu]-1)*(*nVert)+from[segu]-1];
                                   Sigajbj = ILM[(from[segu]-1)*(*nVert)+to[segu]-1];
                                   Sigbjbj = ILM[(to[segu]-1)*(*nVert)+to[segu]-1];
                                   Sigaiai = ILM[(from[segv]-1)*(*nVert)+from[segv]-1];
                                   Sigaibi = ILM[(from[segv]-1)*(*nVert)+to[segv]-1];
                                   Sigbibi = ILM[(to[segv]-1)*(*nVert)+to[segv]-1];
                                   Sigajai = ILM[(from[segu]-1)*(*nVert)+from[segv]-1];
                                   Sigajbi = ILM[(from[segu]-1)*(*nVert)+to[segv]-1];
                                   Sigbjai = ILM[(to[segu]-1)*(*nVert)+from[segv]-1];
                                   Sigbjbi = ILM[(to[segu]-1)*(*nVert)+to[segv]-1];
                       
                                   li = lenEdges[segv];
                                   lj = lenEdges[segu];
                       
                                  vec[id] = (pow(1-tu,2.0) * Sigajaj + pow(tu,2.0) * Sigbjbj + 2*tu*(1-tu) * Sigajbj +
                                             pow(1-tv,2.0) * Sigaiai + pow(tv,2.0) * Sigbibi + 2*tv*(1-tv) * Sigaibi -
                                             2*(1-tu)*(1-tv) * Sigajai - 2*(1-tu)*tv * Sigajbi -
                                             2*tu*(1-tv) * Sigbjai - 2*tu*tv * Sigbjbi +
                                             tv*(1-tv)*li + tu*(1-tu)*lj -
                                             2*aux*fmin(tu*(1-tv),tv*(1-tu))*li);			                                
			         	  
			  }
			  id=id+1;
                            
         }
     }                  
 
}







