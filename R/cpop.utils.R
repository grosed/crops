
########################################################################################
###################### coeff update#####################################################
### NEW VERSION -- UNEVEN LOCATIONS AND VARYING MEAN
###avoids loop
########################################################################################

coeff.update.uneven.var=function(coeffs,S,SXY,SS,SX,SX2,SP,x0,taustar,beta){
  
  coeff.new<-coeffs
  coeff.new[,2]=coeffs[,1] 
  coeff.new[,1]<-taustar
  
  sstar<-coeff.new[,2]
  Xs<-x0[sstar+1]
  Xt<-x0[taustar+1]
  seglen=Xt-Xs





  #A<-(SX2[taustar+1]-SX2[sstar+1]-2*Xs*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xs^2)/(seglen^2)
  #B<- 2*( (Xt+Xs)*(SX[taustar+1]-SX[sstar+1])-(SP[taustar+1]-SP[sstar+1])*Xt*Xs-(SX2[taustar+1]-SX2[sstar+1]))/(seglen^2)
  #C<-(-2)/(seglen)*(SXY[taustar+1]-SXY[sstar+1]-Xs*(S[taustar+1]-S[sstar+1]))
  #D<- (SS[taustar+1]-SS[sstar+1])
  #E<-(-2)/(seglen)*(Xt*(S[taustar+1]-S[sstar+1])-(SXY[taustar+1]-SXY[sstar+1]))
  #FF<-(SX2[taustar+1]-SX2[sstar+1]-2*Xt*(SX[taustar+1]-SX[sstar+1])+(SP[taustar+1]-SP[sstar+1])*Xt^2)/(seglen^2)


  coeffs.cpp<-coeffs_update_cpp(SXY,SX2,S,SS,Xs,SX,SP,seglen,taustar,sstar,Xt)
  A<-coeffs.cpp[[1]]
  B<-coeffs.cpp[[2]]
  C<-coeffs.cpp[[3]]
  D<-coeffs.cpp[[4]]
  E<-coeffs.cpp[[5]]
  FF<-coeffs.cpp[[6]]



  m=length(sstar)
  ind1=(1:m)[FF==0 & coeffs[,3]==0 & B==0]
  ind2=(1:m)[FF==0 & coeffs[,3]==0 & B!=0]
  ind3=(1:m)[!(FF==0 & coeffs[,3]==0)]
  if(length(ind1)>0){
    coeff.new[ind1,5]<-coeffs[ind1,5]+D[ind1]+beta
    coeff.new[ind1,4]<-C[ind1]
    coeff.new[ind1,3]<-A[ind1] 
  }
  
  if(length(ind2)>0){
    coeff.new[ind2,5]<-coeffs[ind2,5]+(-E[ind2]-coeffs[ind2,4])/B[ind2]+beta
    coeff.new[ind2,4]<-0
    coeff.new[ind2,3]<-0
  }
  
  if(length(ind3)>0){
    coeff.new[ind3,5]<-coeffs[ind3,5]+D[ind3]-(coeffs[ind3,4]+E[ind3])^2/(4*(coeffs[ind3,3]+FF[ind3]))+beta
    coeff.new[ind3,4]<-C[ind3]-(coeffs[ind3,4]+E[ind3])*B[ind3]/(2*(coeffs[ind3,3]+FF[ind3]))
    coeff.new[ind3,3]<-A[ind3]-(B[ind3]^2)/(4*(coeffs[ind3,3]+FF[ind3]))
  }
  
  return(coeff.new)
}


############################################################################################################
##first prune of functions
## xx is matrix of functions
############################################################################################################
prune1=function(xx,taustar){
  min.vals<-xx[,5]-(xx[,4]^2)/(2*xx[,3])
  m=min(min.vals[xx[,2]==taustar-1])
  keep<-c(rep(T,length(min.vals[xx[,2]!=taustar-1])),(min.vals[xx[,2]==taustar-1])<=m)
  return(keep) ##keep only those with a smaller minimum.
}


##########################################################################################################
##second pruning
## again x is matrix of quadratics
##version to avoid nested loops
##########################################################################################################

prune2b=function(x){
  Sets<-list()
  n=length(x[,1])
  vec=(1:n)
  
  tcurr= -Inf
  
  whichfun<-which(x[,3]==min(x[,3])) #which element of vec gives min value at -Infinity--smallest theta^2 coeff; then largest theta coeff; then smallest constant
  whichfun<-whichfun[which(x[whichfun,4]==max(x[whichfun,4]))]
  whichfun<-whichfun[which(x[whichfun,5]==min(x[whichfun,5]))]
  whichfun=whichfun[1] #####Introduced this to avoid multiple minimima
  
  Sets[[whichfun]]<-c(tcurr)
  diffcoeffs=matrix(NA,nrow=n,ncol=3)
  intercepts=rep(NA,n)
  disc=rep(NA,n)
  while(length(vec)>1){ #while functions being considered is bigger than 1
    intercepts[1:n]<-NA
    diffcoeffs[1:(length(vec)),]<-t(t(x[vec,3:5])-x[whichfun,3:5]) #difference between coeffs at i and current function
    disc[1:(length(vec))]<-diffcoeffs[1:(length(vec)),2]^2-4*diffcoeffs[1:(length(vec)),1]*diffcoeffs[1:(length(vec)),3] #discriminent of difference quad
    
    ind1=(1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]==0] ##disc>0 for quadratic to cross.
    ind2=(1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]!=0] ##disc>0 for quadratic to cross.
    
    if(length(ind1)>0){
      r1= - diffcoeffs[ind1,3]/diffcoeffs[ind1,2]
      if(sum(r1>tcurr)>0){
        intercepts[ind1[r1>tcurr]]= r1[r1>tcurr]
      }
    }
    if(length(ind2)>0){
      r1=(-diffcoeffs[ind2,2]-sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
      r2=(-diffcoeffs[ind2,2]+sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
      ##only want roots if > tcurr
      if(sum(r1>tcurr)>0){
        intercepts[ind2[r1>tcurr]]=r1[r1>tcurr]
      }
      if(sum(r1<=tcurr & r2>tcurr)>0){
        intercepts[ind2[r1<=tcurr & r2>tcurr]]=r2[r1<=tcurr & r2>tcurr]  
      }
    }
    
    loggy<-!is.na(intercepts)
    loggy[vec==whichfun]<-T
    if(!sum(!is.na(intercepts))==0){ #if at least one intercept value is not na
      tcurr<-min(intercepts,na.rm=T) #change tcurr to first intercept
      whichfunnew<-vec[which(intercepts==tcurr)[1]] #whichfunnew is set as value which first intercept occurs     
      Sets[[whichfun]]<-c(Sets[[whichfun]],tcurr) #add intercept to current function opt interval (to close it)
      if(whichfunnew>length(Sets)){Sets[[whichfunnew]]<-c(tcurr)}else{
        Sets[[whichfunnew]]<-c(Sets[[whichfunnew]],tcurr)} #add intercept to new fucntion interval (opening it)
      whichfun<-whichfunnew #change current function to new function
      
    }
    vec<-vec[loggy[(1:length(vec))]]
    
  }
  Sets[[whichfun]]<-c(Sets[[whichfun]],Inf)
  
  
  output1 <- do.call(rbind,lapply(Sets,length))
  
  return(which(output1[,1]!=0))
}

#############################################################################
#########second pruning in C##################################################
####################################################################################

prune2.c<-function(x){
  nrows<-dim(x)[[1]]
  # coutput<- .C("prune2R",as.double(as.vector((t(x)))),as.integer(nrows),Sets=integer(nrows))
  coutput<-prune(as.vector((t(x))),nrows)
  #result<-which(coutput$Sets==1)
  result<-which(coutput==1)
  return(result)
}

###########################################################################################################
###################################PELT pruning function###################################################
###########################################################################################################

peltprune=function(x,beta){
  if(length(x)==5) return(1)
  minx<-x[,5]
  index=(1:dim(x)[1])[x[,3]>0]
  if(length(index)>0) minx[index]<-x[index,5]-x[index,4]^2/(4*x[index,3]) ####NEED THIS TO BE x[,5] if x[,4]=x[,3]=0
  return(which(minx<=(min(minx)+2*beta)))  
}

###########################################################################################################
##################coverts null to na#######################################################################
###########################################################################################################

null2na<-function(vec){
  if(is.null(vec)){vec<-NA}
  return(vec)
}