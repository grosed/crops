crops.impl <- function(f,beta_star,max_iterations)
{
  n <- 0 
  res <- set()
  while(!set_is_empty(beta_star) && n < max_iterations)
  {
     beta <- as.list(beta_star)[[1]]
     beta_0 <- beta[[1]]
     beta_1 <- beta[[2]]
     Qm_0 <- f(beta_0)[[1]]
     Qm_1 <- f(beta_1)[[1]]
     m_0 <- f(beta_0)[[2]]
     m_1 <- f(beta_1)[[2]]
     if(m_0 > m_1 + 1)
     {
        beta_int <- (Qm_1 - Qm_0)/(m_0-m_1)
	m_int <- f(beta_int)[[2]]
	if(m_int != m_1 && m_int != m_0)
	    {   
		beta_star <- set_union(beta_star,set(tuple(beta_int,beta_1)),set(tuple(beta_0,beta_int)))
	    }
     }
     beta_star <- beta_star - set(beta)
     res <- set_union(res,beta)
     n <- n + 1
  }
  return(res)
}


#' Generic implementation of the crops algorithm (ref goes here).
#'
#' @name crops 
#'
#' @description Provides a generic implementation of the crops (changepoints for a range of penalties) algorithm of Haynes et al. (2014)  which efficiently searches a range of penalty values in multiple changepoint problems.
#' The crops algorithm finds the optimal segmentations for a different number of segments without incurring as large a computational cost as solving the constrained optimisation problem
#' for a range of values for the number of changepoints. To make the method generic, the user must provide a function that maps a penalty value to the results obtained by a penalised cost
#' changepoint method, and formats these results in a specific way. This interface to the generic method is similar to that as used by the \pkg{optimx} package.  
#'
#' @param method A function mapping a penalty value to the results obtained by a penalised cost changepoint method. The function must return a list containing the cost and
#' number of changepoints corresponding to the optimal segmentation as determined by a penalised cost changepoint method. Additional items can be appended to list which get
#' propagated through the call to \code{crops} 
#' @param beta_min A positive numeric value indicating the smallest penalty value to consider.
#' @param beta_max A positive numeric value indicating the maximum penalty value to consider.
#' @param max_iterations Positive non zero integer. Limits the maximum number of iterations of the crops algorithm to \code{max_iterations}. Default value is \code{max_iterations=Inf}.
#' @param ... Additional parameters to pass to the underlying changepoint method if required.
#'
#' @return An instance of an S4 class of type \code{crops.class}. 
#'
#' @references \insertRef{crops-article}{crops}
#' @references \insertRef{optimx-1}{crops}
#' @references \insertRef{optimx-2}{crops}
#' @references \insertRef{optimx-package}{crops}
#' @references \insertRef{fpop-article}{crops}
#' @references \insertRef{fpop-package}{crops}
#'
#' @examples
#' library(cpop)
#' library(crops)
#' 
#' # create a function to wrap a call to cpop for use with crops
#' cpop.crops.eg <- function(y,x,grid,sd,minseglen,prune.approx)
#' {
#'   return(function(beta)
#'          {
#'            res <- cpop(y=y,x=x,grid=grid,beta=beta,sd=sd,minseglen=minseglen,
#'                        prune.approx=prune.approx)
#'            cpts <- changepoints(res)
#'            ncpts <- nrow(cpts)
#'            cost <- cost(res) - ncpts*beta
#'            return(list(cost,ncpts,res))
#'          }
#'        )
#' }
#' 
#' # generate some data
#' set.seed(0)
#' x <- seq(0,1,0.01)
#' n <- length(x)
#' sd <- rep(0.1,n)
#' mu <- c(2*x[1:floor(n/2)],2 - 2*x[(floor(n/2)+1):n])
#' y <- rnorm(n,mu,sd)
#' 
#' # calculate the changepoint locations and cost over a range of penalty values
#' res <- cpop.crops(y,x,sd=sd,beta_min=0.4*log(length(y)),beta_max=2.5*log(length(y)))
#' 
#' # summarise results
#' summary(res)
#' 
#' # extract the segmentations
#' segmentations(res)
#' 
#' # plot the results
#' plot(res)
#'
crops <-
function(method,beta_min,beta_max,max_iterations=Inf,...)
    {
       # appease package checks
       . <- NULL
       check_crops_arguments(method,beta_min,beta_max,max_iterations)
       CPT <- (. %>% method(...) %T>% check_method_return_values %>% as.tuple) %>% memoise
       res <- crops.impl(CPT,set(tuple(beta_min,beta_max)),max_iterations)
       # res <- crops.impl.2(CPT,set(tuple(beta_min,beta_max)),max_iterations)
       object <- crops.class(CPT,res %>% unlist %>% as.set)
       return(object)
    } 




