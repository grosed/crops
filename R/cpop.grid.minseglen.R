

cpop.grid.minseglen<-function(y,x=NULL,grid=NULL,minseg=0,beta,sigsquared=1,printiteration=FALSE,PRUNE.APPROX=FALSE)
{


CPOP.grid.minseg<-function(y,x=NULL,grid=NULL,minseg=0,beta,sigsquared=1,printiteration=FALSE,PRUNE.APPROX=FALSE){

  n<-length(y)
  if(is.null(x)) x= 1:n
  if(is.null(grid)) grid=x  ###SET grid to x if not defined
  if(max(grid)<max(x)) grid=c(grid,max(x))
  if(length(sigsquared)!=n) sigsquared=rep(sigsquared[1],n)
  
  if(minseg<0) minseg=0
  
  ngrid=length(grid) ##ngrid is the number of grid points -- below n becomes ngrid throughout
  
  ##CHANGE HERE -- x0 is grid with an addition of a 0th grid point. Defined in terms of grid not x
  x0=c(2*x[1]-x[2],grid)
  
  if(minseg*2> max(x0)-min(x0)) stop("Minseg is greater than half the grid width => no changepoints possible")
  
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
    if(grid[taustar]-x0[1]>=minseg){ ##ONLY START ADDING CHANGES IF FIRST SEG > minseg
      ##check on which values of most recent change are possible
      index=(1:length(CPvec))[grid[taustar]-x0[coeffs[,1]+1]>=minseg]     
      new.CPvec<-paste(CPvec[index],taustar,sep=",")
      ##update coefficients 
      new.coeffs=coeff.update.uneven.var(coeffs[index,,drop=FALSE],S,SXY,SS,SX,SX2,SP,x0,taustar,beta)
      new.coeffs.p=new.coeffs  
      if(taustar!=ngrid){##n->ngrid #skip pruning on last step
    ###################################################pruning bit##########  
        if(length(new.coeffs[,1])>1){
      ##added###
      ###########
        # keep2=prune2b(new.coeffs.p) ##find set of functions to keep
	keep2=prune2.c(new.coeffs.p)
        new.coeffs.p=new.coeffs.p[keep2,]
        new.CPvec=new.CPvec[keep2]
        }
    ####PELT PRUNE############################
        if(PRUNE.APPROX){
          j=which.max(x0[-1]*(x0[taustar+2]-x0[-1]>=minseg))
          ###IDEA IS WE CAN PELT PRUNE BASED ON COEFFS FOR CHANGE AT j =>
          index=(1:length(CPvec))[grid[j]-x0[coeffs[,1]+1]>=minseg]
          if(length(index)>1){
            coeffs.pelt=coeff.update.uneven.var(coeffs[index,,drop=FALSE],S,SXY,SS,SX,SX2,SP,x0,j,beta)
            keeppelt<-peltprune(coeffs.pelt,1.5*beta) ##this holds for new.coeffs by choice above.
            if(length(keeppelt)<length(index)){
              coeffs<-coeffs[-index[-keeppelt],]
              CPvec<-CPvec[-index[-keeppelt]]
            }
          }
         }
        }
    ##########################################
    CPvec<-c(CPvec,new.CPvec) 
    coeffs<-rbind(coeffs,new.coeffs.p)
    }
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


out<-CPOP.grid.minseg(y,x,grid,minseg,beta,sigsquared,FALSE,FALSE)
out$changepoints<-c(x[1],out$changepoints)
return(out)


}