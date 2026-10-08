crops.impl <- function(f,beta_star)
{
  res <- set()
  while(!set_is_empty(beta_star))
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
#' library(crops)
#' 
#' # create a segmentation function for demonstration
#' toy_segmentation_function <- function(beta)
#' {
#'    x <- c(0.0,1.0,1.4,1.5,1.9)
#'    C <- c(2.5,4.0,4.4,5.6,9.5)
#'    n <- c(10,8,8,7,6)
#'    i <- findInterval(beta,x)
#'    return(as.list(c(C[i],n[i])))
#' }
#'
#' # plot the segmentation function
#' beta <- seq(0.0,3.0,0.01)
#' values <- Map(toy_segmentation_function,beta)
#' C <- as.numeric(Map(function(.) (.)[[1]],values))
#' n <- as.numeric(Map(function(.) (.)[[2]],values))
#' plot(beta,n,type="l",col="red",ylim=c(2.0,10))
#' lines(beta,C,type="l",col="green")
#' legend("bottomright",
#'        legend = c("Unpenalized Cost", "Number of Changepoints"),
#'        col = c("green", "red"),
#'        pch = c(19, 17),
#'        lty = c(1, 1))
#' 
#' # analyse with crops
#' results <- crops(toy_segmentation_function,0.0,3.0)
#'
#' # x summarise results
#' Map(results@method,as.list(results@betas))
#'
crops <-
function(method,beta_min,beta_max,...)
    {
       # appease package checks
       . <- NULL
       check_crops_arguments(method,beta_min,beta_max)
       CPT <- (. %>% method(...) %T>% check_method_return_values %>% as.tuple) %>% memoise
       res <- crops.impl(CPT,set(tuple(beta_min,beta_max)))
       object <- crops.class(CPT,res %>% unlist %>% as.set)
       return(object)
    } 




