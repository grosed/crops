
#' Generic implementation of the crops algorithm (ref goes here).
#'
#' @name crops 
#'
#' @description Provides a generic implementation of the crops (changepoints for a range of penalties) algorithm of Haynes et al. (2014)  which efficiently searches a range of penalty values in multiple changepoint problems.
#' The crops algorithm finds the optimal segmentations for a different number of segments without incurring as large a computational cost as solving the constrained optimisation problem
#' for a range of values for the number of changepoints. To make the method generic, the user must provide a function that maps a penalty value to the results obtained by a penalised cost
#' changepoint method, and formats these results in a specific way. This interface to the generic method is similar to that as used by the \pkg{optimx} package.  
#'
#' @param method A function mapping a penalty value to the results obtained by a penalised cost changepoint method. The function must return a list containing the cost, number of changepoints, and
#' a vector of changepoint locations corresponding to the optimal segmentation as determined by a penalised cost changepoint method.
#'
#' @param beta_min A positive numeric value indicating the smallest penalty value to consider.
#' @param beta_max A positive numeric value indicating the maximum penalty value to consider.
#' @param max_iterations Positive non zero integer. Limits the maximum number of iterations of the crops algorithm to \code{max_iterations}. Default value is \code{max_iterations=Inf}
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
#' # generate some simple data
#' set.seed(1)
#' N <- 100
#' data.vec <- c(rnorm(N), rnorm(N, 2), rnorm(N))
#'
#' # example one - calling fpop via crops using global scope
#' # need the fpop library
#' library(pacman)
#' p_load(fpop)
#' # create a function to wrap a call to fpop for use with crops
#' fpop.for.crops<-function(beta)
#'     {
#'        # Note - this code is taken from the example in the fpop package
#'        fit <- Fpop(data.vec, beta)
#'        end.vec <- fit$t.est
#'        change.vec <- end.vec[-length(end.vec)]
#'        start.vec <- c(1, change.vec+1)
#'        segs.list <- list()
#'        for(seg.i in seq_along(start.vec))
#'            {
#'             start <- start.vec[seg.i]
#'             end <- end.vec[seg.i]
#'             seg.data <- data.vec[start:end]
#'             seg.mean <- mean(seg.data)
#'             segs.list[[seg.i]] <- data.frame(
#'                                     start, end,
#'                                     mean=seg.mean,
#'                                     seg.cost=sum((seg.data-seg.mean)^2))
#'             }
#'         segs <- do.call(rbind, segs.list)
#'         return(list(sum(segs$seg.cost),nrow(segs)-1,segs$end[-length(segs$end)]))
#'     }
#'
#' # now use this wrapper function with crops
#' res<-crops(fpop.for.crops,0.5*log(300),2.5*log(300))
#' # print summary of analysis
#' summary(res)
#' # summarise the segmentations
#' segmentations(res)
#' # visualise the segmentations
#' plot(res)
#' # overlay the data on the segmentations
#' df <- data.frame("x"=1:300,"y"=data.vec)
#' plot(res,df) 
#'
crops <-
function(method,beta_min,beta_max,max_iterations=Inf,...)
    {
       check_crops_arguments(method,beta_min,beta_max,max_iterations)
       CPT<-memoise(function(.)
       {
	   res <- method(.,...) 
           check_method_return_values(res)
           return(as.tuple(res))
       })
       tryCatch(
       {
       beta_star <- set(tuple(beta_min,beta_max))
       log <- set()
       iterations <- 0
       while(!set_is_empty(beta_star) && iterations < max_iterations)
           {
             for(t in beta_star)
                 {
                   beta_0 <- t[[1]]
                   beta_1 <- t[[2]]
                   left <- as.tuple(c(beta_0,CPT(beta_0)))
                   right <- as.tuple(c(beta_1,CPT(beta_1)))
                   log <- set_union(log,set(left),set(right))
                   if(left[[3]] > right[[3]] + 1)
                    {
                        beta_int <- (right[[2]] - left[[2]]) / (left[[3]] - right[[3]])
                        int <- as.tuple(c(beta_int,CPT(beta_int)))
                        log <- set_union(log,set(int))
                        if(int[[3]] > right[[3]])
                        {
                           beta_star <- set_union(beta_star,set(tuple(beta_0,beta_int)),set(tuple(beta_int,beta_1)))
                        }
                    }
                  beta_star <- set_symdiff(beta_star,set(t))
                 }
	     iterations <- iterations + 1	 
            }
	if(iterations == max_iterations)
	{
	   warning(paste("maximum number of iterations (=",max_iterations,") reached in crops.",'\n',sep=""))
	}
        return(crops.class(log,beta_min,beta_max,iterations))
	},
	interrupt = function(e)
	            {
		      warning("crops interrupted via CTRL-C. Expect results to be incomplete and/or corrupted.")
		      return(crops.class(log,beta_min,beta_max,iterations))
	            }
	)
    } 




	  