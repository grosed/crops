setOldClass("memoised")
setOldClass("set")
.crops.class<-setClass("crops.class",representation(method="memoised",betas="set"))
crops.class<-function(method,betas)
{
    .crops.class(method=method,betas=betas)
}


#' Summary of segmentations by penalty value
#'
#' @name segmentations
#'
#' @description Produces a summary of the segmentations for each penalty value in the form of a data frame.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{crops}}.
#' 
#' @return A data frame containing the penalties, costs, penalised costs, and changepoint locations
#'
#' @rdname segmentations-methods
#'
#' @aliases segmentations,crops.class-method
#' 
#' @seealso \code{\link{crops}}.
#'
#' @examples
#' # see the crops example
setGeneric("segmentations",function(object) {standardGeneric("segmentations")})
setMethod("segmentations",signature=list("crops.class"),
          function(object)
          {
	    # appease package checks
            . <- NULL
            segs <- Map(object@method,unlist(object@betas))
            n <- segs %>% Map(function(.) .[[2]],.) %>% Map(length,.) %>% unlist %>% max
            mat <- segs %>% 
                   Map(function(.) .[[2]],.) %>% 
                   Map(function(.) c(.,rep(NA,n-length(.))),.) %>%
                   Reduce(rbind,.,matrix(nrow=0,ncol=n),right=TRUE)
                   colnames(mat) <- Map(function(.) paste("cpt.",.,sep=""),1:n) %>% unlist      
            return(      
                    Map(function(beta,seg) tibble(beta=beta,Qm=seg[[1]],Q=seg[[1]]+beta*length(seg[[2]]),
                        m=length(seg[[2]])),
                        unlist(object@betas),
                        segs) %>%                                          
                    Reduce(add_row,.,tibble(beta=numeric(),Qm=numeric(),Q=numeric(),m=numeric())) %>%
                    cbind(.,mat)
                   )           
})


#' Pretty printing for crops results
#'
#' @name print
#'
#' @description Pretty prints a summary of a crops result
#'
#' @docType methods
#'
#' @param  x An instance of an S4 class produced by \code{\link{crops}}.
#'
#' @rdname print-methods
#'
#' @aliases print,crops.class-method
#' 
#' @seealso \code{\link{crops}}.
#'
#' @examples
#' # see the crops example
#'
setMethod("print",signature=list("crops.class"),
          function(x)
          {
	    summary(x)            
})


#' Visualisation of data, costs, penalty values and changepoint locations.
#'
#' @name plot
#'
#' @description Plot methods for an S4 object returned by \code{\link{crops}}. The plot can also be combined with the original data if required.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{crops}}.
#' @param y A dataframe containing the locations and values of the data points. The data plot is plotted below, and is aligned with, the changepoint locations.
#' 
#' @return A ggplot object.
#'
#' @rdname plot-methods
#'
#' @aliases plot,crops.class,data.frame-method
#' 
#' @seealso \code{\link{crops}}.
#'
#' @examples
#' # see the crops example
#'
#' @export  
setMethod("plot",signature=list("crops.class","data.frame"),
          function(x,y)
          {
	    p1 <- plot(x)
	    p2 <- ggplot(data=y) + geom_line(aes(x=x,y=y))
            return(plot_grid(p1,p2,ncol=1))
          })

#' @name plot
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{crops}}.
#' 
#' @rdname plot-methods
#'
#' @aliases plot,crops.class,missing-method
#'
#' @examples
#' # see the crops example
#'
#' @export
setMethod("plot",signature=list("crops.class","missing"),
          function(x)
          {
	     # appease ggplot and tidyverse
	     . <- Q <- Qm <- m <- value <- dummy <- NULL
	     object <- x
             df <- segmentations(object)
	     df <- cbind(df,data.frame("dummy"=1:nrow(df)))
             p <- df %>%
	          subset(.,select = -c(beta,Q,Qm,m)) %>%
                  melt(., id=c("dummy")) %>% 
                  .[complete.cases(.), ] %>%
                  ggplot(.,aes(x=value,y=dummy)) %>% 
                  add(geom_point()) %>%
                  add(labs(x="location",y="penalty")) %>%
                  add(geom_hline(aes(yintercept=dummy))) %>%
                  add(scale_y_continuous(breaks = seq(1:nrow(df)),labels=signif(df$beta,digits=3),sec.axis = sec_axis( ~.,breaks = seq(1:nrow(df)),labels=signif(df$Q,digits=4),name="penalised cost"))) %>%
                  add(theme_bw())
             return(p)       
          })


#' Summary of crops result
#'
#' @name summary
#'
#' @description Prints a short summary of a crops result.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{crops}}.
#'
#' @rdname summary-methods
#'
#' @aliases summary,crops.class-method
#' 
#' @seealso \code{\link{crops}}.
#'
#' @examples
#' # see the crops example
#'
setMethod("summary",signature=list("crops.class"),
          function(object)
          {
            cat("crops analysis",sep="")
	    cat('\n',sep="")
	    cat('\n',sep="")
	    cat("minimum penalty value = ",min(object@betas)," : maximum penalty value = ",max(object@betas),sep="")
	    cat('\n',sep="")
	    segs <- segmentations(object)
	    cat("number of segmentations calculated : ",nrow(segs),sep="")
	    cat('\n',sep="")	    
	    cat("least number of changepoints  = ",min(segs$m), " : maximum number of changepoints = ",max(segs$m),sep="")
	    cat('\n',sep="")
        invisible()
})


#' Remove duplicate entries from a crops result
#'
#' @name unique
#'
#' @description Removes duplicate entries from a crops result. A duplicate entry is one having the same number of changepoints as another entry.
#' Note that the changepoint locations and the associated penalty and cost values are not taken into consideration. The \code{unique} function can be useful
#' for simplifying plots and the details produced by \code{segmentations}.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{crops}}.
#'
#' @return An instance of the S4 class type \code{crops.class}. This is the same type as produced by the \code{\link{crops}} function.
#'
#' @rdname unique-methods
#'
#' @aliases unique,crops.class-method
#' 
#'
setMethod("unique",signature=list("crops.class"),
         function(x)
             {
	     	# appease package checks
                . <- NULL
                object<-x
                hash_map <- new.env()
                keys <- object@betas %>% 
                unlist %>% 
                Map(object@method,.) %>% 
                Map(function(.) .[2],.) %>% 
                Map(as.character,.)
                key_value_pairs <- Map(tuple,keys,object@betas %>% unlist)
                hash_map <-    
                key_value_pairs %>%  
                Reduce(function(pair,map) {map[[pair[[1]]]] <- pair[[2]]; return(map);},
                   .,
                   hash_map,right=TRUE)
                object@betas <-    
                Map(function(key) hash_map[[key]],
                    hash_map %>% ls) %>%
                unname %>% 
                as.set
                return(object)
             })


#' Subset crops results based on penalty values
#'
#' @name subset
#'
#' @description Removes entries from a crops result that fall outside a specified range of penalty values. 
#' The \code{subset} function can be useful for simplifying plots and the details produced by \code{segmentations}.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{crops}}.
#' @param beta_min A positive numeric value specifying the minimum penalty value for entries in the crops result. Default value is 0.
#' @param beta_max A positive numeric value specifying the maximum penalty value for entries in the crops result. Default value is Inf.
#'
#' @return An instance of the S4 class type \code{crops.class}. This is the same type as produced by the \code{\link{crops}} function.
#'
#' @rdname subset-methods
#'
#' @aliases subset,crops.class-method
#' 
#'
setMethod("subset",signature=list("crops.class"),
         function(x,beta_min=0,beta_max=Inf)
             {
	     	# appease package checks
                . <- NULL
	        object <- x
                object@betas %<>% 
                unlist %>% 
                Filter(function(.) . <= beta_max & . >= beta_min,.) %>% 
                as.set
                return(object)            
             })


