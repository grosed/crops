#' A technique for ...
#'
#' @name crops 
#'
#' @description A technique for ... 
#'
#' @param method A ...
#' @param beta_min A ...
#' @param beta_max A ...
#' @param max_iterations A ...
#' 
#' @return An instance of an S4 class of type crops.class. 
#'
#' @examples
#' # generate some simple data
#' set.seed(1)
#' N <- 100
#' data.vec <- c(rnorm(N), rnorm(N, 2), rnorm(N))
#'
#' # example one - calling fpop with crops using global scope
#' # need the fpop library
#' library(pacman)
#' p_load(fpop)
#' # create a function to wrap a call to fpop for use with crops
#' fpop.for.crops<-function(beta)
#'     {
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
#'         return(list(sum(segs$seg.cost),nrow(segs),segs$end))
#'     }
#'
#' # now use this wrapperfunction with crops
#' res<-crops(fpop.for.crops,0.5*log(300),2.5*log(300))
#' # and plot the results
#' plot(res)
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
        return(crops.class(log))
	},
	interrupt = function(e)
	            {
		      warning("crops interrupted via CTRL-C. Expect results to be incomplete and/or corrupted.")
		      return(crops.class(log))
	            }
	)
    } 




	  