empVar <- function(nlags,hmax,tol,zVec,distMat){

 h = seq(0,hmax,l=nlags);
 vari=c();
 madog=c();
 for(i in 1:nlags){
    ind = which(distMat<h[i]+tol & distMat>h[i]-tol,arr.ind=TRUE);
    ind = ind[ind[,2]>=ind[,1],];
    aux=zVec[ind[,1]]-zVec[ind[,2]];
    vari[i] = mean(aux^2);
    madog[i] = mean(abs(aux));
 }
 return(rbind(vari,madog))
}
