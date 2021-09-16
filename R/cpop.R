
.cpop.class<-setClass("cpop.class",representation(changepoints="numeric",
						  y="numeric",
						  x="numeric",
						  beta="numeric",
						  sd="numeric"))

cpop.class<-function(y,x,beta,sd,changepoints)
{
    .cpop.class(y=y,
                x=x,
		changepoints=changepoints,
		beta=beta,
		sd=sd)
}



#' A function for calculating the cost of a model fitted by cpop
#'
#' @name cost
#'
#' @description Calculates the cost of a model fitted by cpop using the residual sum of squares and the penalty values
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @rdname cost-methods
#'
#' @aliases cost,cpop.class-method
#'
#' @export
setGeneric("cost",function(object) {standardGeneric("cost")})
setMethod("cost",signature=list("cpop.class"),
          function(object)
          {
            df <- fitted(object)
	    return(sum(residuals(object)^2/object@sd^2)+object@beta*(length(object@changepoints)-2))
          })	      


#' A function for generating simulated data
#'
#' @name simulate
#'
#' @description Generates simulated data for use with \code{\link{cpop}}.
#' 
#' @param x A numeric vector of containing the locations of the data.
#' @param changepoints A numeric vector of changepoint locations.
#' @param change.slope A numeric vector indicating the change in slope at each changepoint. The initial slope is assumed to be 0.
#' @param sigma The residual standard deviation. Can be a single numerical value or a vector of values for the case of varying residual standard deviation.Default value is 1.
#' @return A vector of simulated y values.
#'
#' @examples
#' library(cpop)
#' set.seed(1)
#' changepoints=c(0,25,50,100)
#' change.slope=c(0.2,-0.3,0.2,-0.1)
#' x=1:200
#' sig=1+x/200
#' y<-simulate(x,changepoints,change.slope,sig)
#' res<-cpop(y,x,beta=2*log(length(x)),sd=sig)
#' summary(res)
#' plot(res)
#'
#' @export
simulate<-function(x,changepoints,change.slope,sigma=1)
{
  K=length(changepoints)
  mu=rep(0,length(x))
  for(k in 1:K)
  {
     mu=mu+change.slope[k]*pmax(x-changepoints[k],0)
  }
  y=mu+rnorm(length(x),0.0,sigma)
  return(y)
}


#' Visualisation of Changepoint Locations and Data
#'
#' @name plot
#'
#' @description Plot methods for S4 objects returned by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param x An instance of an cpop S4 class produced by \code{\link{cpop}}.
#' 
#' @return A ggplot object.
#'
#' @rdname plot-methods
#'
#' @aliases plot,cpop.class-method
#' 
#' @export 
setMethod("plot",signature=list("cpop.class"),function(x)
{
  object <- x
  df <- data.frame("x"=object@x,"y"=object@y)
  cpts<-object@changepoints
  # appease ggplot2
  y <- NULL
  p <- ggplot(data=df, aes(x=x, y=y))
  p <- p + geom_point(alpha=0.3)
  if(length(cpts) > 2)
  {
   for(cpt in cpts[2:(length(cpts)-1)])
   {
     p <- p + geom_vline(xintercept = cpt,color="red")
   }
  }
  # appease ggplot2
  x0 <- y0 <- x1 <- y1 <- NULL  
  #df<-data.frame("xs"=obj@x[cpts][1:(length(obj@x[cpts])-1)],
  #	         "ys"=obj@y_hat[cpts][1:(length(obj@y_hat[cpts])-1)],
  #		 "xends"=obj@x[cpts][2:length(obj@x[cpts])],
  #		 "yends"=obj@y_hat[cpts][2:length(obj@y_hat[cpts])])
  p <- p + geom_segment(data=fitted(object),aes(x=x0,y=y0,xend=x1,yend=y1))
  p <- p + theme_bw()
  return(p)
})



#' Summary of cpop Analysis.
#'
#' @name summary
#'
#' @description Summary method for results produced by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @rdname summary-methods
#'
#' @aliases summary,cpop.class-method
#'
#' @export
setMethod("summary",signature=list("cpop.class"),function(object)
{
  cat('\n',"cpop analysis with n = ",length(object@x)," and penalty (beta)  = ",object@beta,'\n\n',sep="")
  if(length(object@changepoints) == 2)
  {
    cat("No changepoints detected",'\n\n',sep="")
  }
  else
  {
    msg<-paste(length(object@changepoints)-2," ",sep="")
    if(length(object@changepoints) == 3)
     {
       msg<-paste(msg," changepoint",sep="")
     }
     else
     {
       msg<-paste(msg," changepoints",sep="")
     }
     msg<-paste(msg," detected at x = ",'\n',sep="")
     cat(msg)
     msg<-""
     for(cpt in object@changepoints[2:(length(object@changepoints)-1)])
     {
       msg<-paste(msg,cpt,sep=" ")
     }
     msg<-paste(msg,'\n',sep="")
     cat(msg)
  }
  df <- fitted(object)
  cat("fitted values : ",'\n',sep="")
  print(df)
  cat('\n',"overall RSS = ",sum(df$RSS),'\n',sep="")
  cat("cost = ",cost(object),'\n',sep="")
  invisible()
})


#' Displays an S4 object produced by cpop.
#'
#' @name show
#'
#' @description Displays an S4 object produced by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @rdname show-methods
#'
#' @aliases show,cpop.class-method
#' 
#'
#' @export
setGeneric("show",function(object) {standardGeneric("show")})
setMethod("show",signature=list("cpop.class"),function(object)
{
    summary(object)
    invisible()
})


#' Extract Model Fitted Values
#'
#' @name fitted
#'
#' @description Extracts the fitted valus produced by \code{\link{cpop}}.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#'
#' @return A data frame containing the endpoint coordinates for each line segment fitted between
#' the detected changepoints. The data frame also contains the gradient and intercept values
#' for each segment and the corresponding residual sum of squares (RSS)..
#'
#' @rdname fitted-methods
#'
#' @aliases fitted,cpop.class-method
#'
#' @export
setGeneric("fitted",function(object) {standardGeneric("fitted")})
setMethod("fitted",signature=list("cpop.class"),
          function(object)
          {
	  x<-sort(unique(c(object@x,object@changepoints)))
	  y_hat<-estimate(object,x)$y_hat
	  cpts<-object@changepoints
	  y_0<-unlist(Map(function(val) y_hat[which(x==val)],cpts))
  	  df<-data.frame("x0"=cpts[1:(length(cpts)-1)],
			 "y0"=y_0[1:(length(cpts)-1)],
		         "x1"=cpts[2:length(cpts)],
		         "y1"=y_0[2:length(y_0)])
          df <- cbind(df,data.frame("gradient"=(df$y1 - df$y0)/(df$x1 - df$x0)))
	  df <- cbind(df,data.frame("intercept"=df$y1 - df$gradient * df$x1))
	  residuals<-residuals(object)
          rss<-unlist(Map(function(a,b) sum((residuals*residuals)[which(object@x >= a & object@x < b)]),cpts[1:(length(cpts)-1)],cpts[2:length(cpts)]))
	  rss[length(rss)]<-rss[length(rss)]+(residuals*residuals)[length(residuals)]
	  df <- cbind(df,data.frame("RSS"=rss))
	  return(df)
          })	      


#' Changepoint Locations
#'
#' @name changepoints
#'
#' @description Creates a data frame containing the locations of the changepoints in terms of the index of the data and the value of the location at that index.
#'
#' @docType methods
#'
#' @rdname changepoints-methods
#'
#' @param object  An instance of an cpop S4 class produced by \code{\link{cpop}}.
#' 
#' @return A data frame.
#' 
#' @aliases changepoints,cpop.class-method
#'
#' @examples
#' library(cpop)
#' # generate some test data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sigma <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sigma)
#'
#' # use the locations in x
#' res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
#' changepoints(res)
#'
#' @export
if(!isGeneric("changepoints")) {setGeneric("changepoints",function(object) {standardGeneric("changepoints")})}
setGeneric("changepoints",function(object) {standardGeneric("changepoints")})
setMethod("changepoints",signature=list("cpop.class"),
          function(object)
          {
	      if(length(object@changepoints) > 2)
	      {
	      	      df <- data.frame("location"=object@changepoints[2:(length(object@changepoints)-1)])	
	      }
	      else
	      {
	      	      df <- data.frame("location"=numeric(0))	
	      }
	      return(df)
          })	      


#' cpop
#'
#'  Algorithm for finding the best segmentation of data for a change-in-slope model.
#' 
#' @param y A vector of length n containing the data.
#' @param x A vector of length n containing the locations of y. Default value is NULL, in which case the locations \code{x = 1:length(y)} are assumed.
#' @param grid An ordered vector of possible locations for the change points. If this is NULL, then this is set to x, the vector of times/locations of the data points.
#' @param minseglen The minimum allowable segment length, i.e. distance between successive changepoints. Default is 0.
#' @param prune.approx Only relevant if a minimum segment length is set. If True, cpop will use an approximate pruning algorithm that will speed up computation but may
#' occasionally lead to a sub-optimal solution in terms of the estimate change point locations. If the minimum segment length is 0, then an exact pruning algorithm is possible and is used.
#' @param beta A positive real value for the penalty incurred for adding a changepoint (prevents over-fitting).
#' @param sd Estimate of residual standard deviation. Can be a single numerical value or a vector of values for the case of varying standard deviation. Default value is 1. 
#'
#' @return An instance of an S4 class of type cpop.class.
#'
#' @references \insertRef{doi:10.1080/10618600.2018.1512868}{cpop}
#'
#' @examples
#'
#' library(cpop)
#' # generate some test data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sigma <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sigma)
#'
#' # use the locations in x
#' res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
#' plot(res)
#'
#' # without locations (note explicit paramater names)
#' res <- cpop(y,beta=2*log(length(y)),sd=sigma)
#' plot(res)
#'
#' # stretch the end of the data
#' x[75:101] <- x[75:101] + seq(from=0,by=0.2,length.out=27)
#' res <- cpop(y,x,beta=2*log(length(y)),sd=sigma)
#' plot(res)
#'  
#' @export
cpop<-function(y,x=1:length(y),grid=x,beta=2*log(length(y)),sd=1,minseglen=0,prune.approx=FALSE)
{
    if(length(sd)!=length(y))
    {
      sd=rep(sd[1],length(y))
    }
    sigsquared<-sd^2
    if(minseglen != 0)
    {
      res<-cpop.grid.minseglen(y,x,grid,beta,sigsquared,minseg=minseglen,FALSE,prune.approx)
    }
    else
    {
      res<-cpop.grid(y,x,grid,beta,sigsquared)
    }
    return(cpop.class(y,x,beta,sd,res$changepoints))
}


design<-function(object,x=object@x)
{
  n=length(x)
  cpts<-object@changepoints[-1]
  p=length(cpts)
  X=matrix(NA,nrow=n,ncol=p+1)
  X[,1]=1
  X[,2]=x-x[1]
  if(p>1)
  {
    for(i in 1:(p-1))
    {
      X[,i+2]=pmax(rep(0,n),x-cpts[i])
    }
  }
  return(X)
}

parameters<-function(object)
{
  n=length(object@y)
  cpts<-object@changepoints[-1]	
  p<-length(cpts)
  W<-diag(object@sd^-2)
  X<-design(object)
  XTX<-t(X)%*%W%*%X
  pars<-as.vector(solve(XTX)%*%t(X)%*%W%*%object@y)
  return(pars)
}

#' A function for estimating the fit of a cpop model
#'
#' @name estimate
#'
#' @description Estimates the fit of a cpop model at the specified locations
#'
#' @param object An instance of an S4 class produced by \code{\link{cpop}}.
#' @param x Locations at which the fit is to be estimated. Default value is the x locations at which the cpop object was defined.
#' @param ... Additional arguments.
#'
#' @rdname estimate-methods
#'
#' @aliases estimate,cpop.class-method
#'
#' @export
setGeneric("estimate",function(object,x,...) {standardGeneric("estimate")})
setMethod("estimate",signature=list("cpop.class"),
          function(object,x=object@x)
          {
             return(data.frame("x"=x,"y_hat"=design(object,x)%*%parameters(object)))
          })	      



residuals<-function(object)
{
  object@y-design(object)%*%parameters(object)
}


cpop.fit<-function(y,x,out.changepoints,sigsquared)
{
  n=length(y)
  if(length(sigsquared)!=n)
  {
     sigsquared=rep(sigsquared[1],n)
  }
  p=length(out.changepoints)
  W=diag(sigsquared^-1)
  X=matrix(NA,nrow=n,ncol=p+1)
  X[,1]=1
  X[,2]=x-x[1]
  if(p>1)
  {
    for(i in 1:(p-1))
    {
      X[,i+2]=pmax(rep(0,n),x-out.changepoints[i])

    }
  }
  XTX=t(X)%*%W%*%X
  beta=as.vector(solve(XTX)%*%t(X)%*%W%*%y)
  fit=X%*%beta
  residuals=y-fit
  return(list(fit=fit,residuals=residuals,X=X,pars=beta))
}








