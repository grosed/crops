cpop.grid<-function(y,x,grid,beta,sigsquared)
{


CPOP.grid.impl<-function(y,x,grid=NULL,beta,sigsquared=1,printiteration=FALSE){

  n<-length(y)
  if(is.null(x)) x= 1:n
  if(is.null(grid)) grid=x  ###SET grid to x if not defined
  if(max(grid)<max(x)) grid=c(grid,max(x))
  if(length(sigsquared)!=n) sigsquared=rep(sigsquared[1],n)
  
  ngrid=length(grid) ##ngrid is the number of grid points -- below n becomes ngrid throughout
  
  ##CHANGE HERE -- x0 is grid with an addition of a 0th grid point. Defined in terms of grid not x
  x0=c(2*x[1]-x[2],grid)
  
  ##THESE UPDATES ARE CHANGED.
  S<-0
  SS<-0
  SX<-0
  SX2<-0
  SXY<-0
  SP<-0
  for(i in 1:ngrid){
    ##set index to be which observations are in the interval
    index=(1:n)[x>x0[i]&x<=x0[i+1]]
    if(length(index)>0){##change in summaries -- sum over observations
      S[i+1]<-S[i]+sum(y[index]/sigsquared[index])
      SS[i+1]<-SS[i]+sum(y[index]^2/sigsquared[index])
      SX[i+1]<-SX[i]+sum(x[index]/sigsquared[index])
      SX2[i+1]<-SX2[i]+sum(x[index]^2/sigsquared[index])
      SXY[i+1]<-SXY[i]+sum(x[index]*y[index]/sigsquared[index])
      SP[i+1]<-SP[i]+sum(1/sigsquared[index])
    }else{ ##no observations so no change in summaries
      S[i+1]<-S[i]
      SS[i+1]<-SS[i]
      SX[i+1]<-SX[i]
      SX2[i+1]<-SX2[i]
      SXY[i+1]<-SXY[i]
      SP[i+1]<-SP[i]
    }
  }
  #######
  
  coeffs<-matrix(0,ncol=5,nrow=1) #first two columns are current time point and most recent changepoint, final three are coefficients for cost
  coeffs[1,5]<--beta
  coeffs[1,1:2]<-c(0,0)
  CPvec<-c("0") #vector storing changepoint values, not used in code but required as an output 
  
  for(taustar in 1:ngrid){ ##n->ngrid is the number of iterations
    new.CPvec<-paste(CPvec,taustar,sep=",")
    ##update coefficients --THIS HAS BEEN CHANGED FROM CPOP CODE
    ##CURRENTLY CHANGE ONLY IN R CODE VERSION
    new.coeffs=coeff.update.uneven.var(coeffs,S,SXY,SS,SX,SX2,SP,x0,taustar,beta)
    new.coeffs.p=new.coeffs
    if(taustar!=ngrid){##n->ngrid #skip pruning on last step
    ###################################################pruning bit##########  
      if(length(new.coeffs[,1])>1){
      ##added###
      ###########
       keep2=prune2.c(new.coeffs.p)
       # keep2=prune2b(new.coeffs.p) ##find set of functions to keep
       new.coeffs.p=new.coeffs.p[keep2,]
       new.CPvec=new.CPvec[keep2]
      }
    ####PELT PRUNE############################
      if(length(new.coeffs)>5 ){
        keeppelt<-peltprune(new.coeffs,beta)
        coeffs<-coeffs[keeppelt,]
        CPvec<-CPvec[keeppelt]
      }
    }
    ##########################################
    CPvec<-c(CPvec,new.CPvec) #prunes both CPvec vector and coeffs matrix
    coeffs<-rbind(coeffs,new.coeffs.p)
     #####################################################
    if(printiteration==TRUE){
      if(taustar%%100==0) cat("Iteration ",taustar,"\n")}
      else if(printiteration!=FALSE){stop("printiteration must be a TRUE or FALSE value")}
    }
  coeffscurr<-coeffs[coeffs[,1]==ngrid,] #matrix of coeffs for end time t=n ##n->ngrid
  if(!is.matrix(coeffscurr)){coeffscurr<-t(as.matrix(coeffscurr))} #makes sure coeffscurr is in the right format
  ##HACK here to make ttemp correct if coeffscurr[,3]==0
  ttemp<-coeffscurr[,5]
  index=(1: (length(coeffscurr)/5))[coeffscurr[,3]>0] 
  if(length(index)>0){ 
    ttemp[index]<-coeffscurr[index,5]-(coeffscurr[index,4]^2)/(4*coeffscurr[index,3])}
  mttemp<-min(ttemp)
  num<-which(ttemp==mttemp)
  
  CPveccurr<-CPvec[coeffs[,1]==ngrid] ##n->ngrid
  CPS<-eval(parse(text=paste("c(",CPveccurr[num],")")))

  return(list(min.cost=mttemp,changepoints.index=CPS,changepoints=grid[CPS])) #return min cost and changepoints
}


   out<-CPOP.grid.impl(y,x,grid,beta,sigsquared,FALSE)
   out$changepoints<-c(x[1],out$changepoints)
   return(out)

}