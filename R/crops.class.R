setOldClass("set")
.crops.class<-setClass("crops.class",representation(log="set",beta_min="numeric",beta_max="numeric",iterations="numeric"))
crops.class<-function(log,beta_min,beta_max,iterations)
{
    .crops.class(log=log,beta_min=beta_min,beta_max=beta_max,iterations=iterations)
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
#' @return A data frame contianing the penalties, costs, penalised costs, and changepoint locations
#'
#' @rdname segmentations-methods
#'
#' @aliases segmentations,crops.class-method
#' 
#' @seealso \code{\link{crops}}.
#'
#' @examples
#' # see the crops example
#'
setGeneric("segmentations",function(object) {standardGeneric("segmentations")})
setMethod("segmentations",signature=list("crops.class"),
          function(object)
          {
            n <- max(unlist(Map(function(.) length(unlist(.)),as.list(object@log)))) - 3
            cpts <- Map(function(y) c(y[[4]],rep(NA,n-length(y[[4]]))), as.list(object@log))
            mat <- Reduce(rbind,cpts,matrix(nrow=0,ncol=n),right=TRUE)
            colnames(mat) <- unlist(Map(function(i) paste("cpt.",i,sep=""),1:n))
            df<-data.frame(beta=numeric(),Qm=numeric(),Q=numeric(),m=numeric())
            for(item in object@log)
               {
                  df[nrow(df)+1,] <- c(item[1:2],item[2] + item[3]*item[1],item[3])
               }
            df <- cbind(df,mat)
            return(df)
            
})


#' Visualisation of data, costs, penatly values and changepoint locations.
#'
#' @name plot
#'
#' @description Plot methods for an S4 object returned by \code{\link{crops}}. The plot can also overlay the original data over the changepoint locaations if required.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{crops}}.
#' @param y A dataframe containing the locations and values of the data points. The data will be overlayed on the changepoint locations.
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
	    # appease ggplot2
	    a <- b <- NULL
	    object <- x
	    data <- y
            p <- plot(object)
            df <- segmentations(object)
            n <- nrow(df)
            y <- data[,2]
            y <- y - min(y)
            y <- 0.75*n*y/(max(y) - min(y))
            y <- y + n/2 - mean(y) 
            x <- data[,1]
            p <- p + geom_line(data=data.frame(a=x,b=y),aes(x=a,y=b),alpha=0.2)
            return(p)
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
	    cat("minimum penalty value = ",object@beta_min," : maximum penalty value = ",object@beta_max,sep="")
	    cat('\n',sep="")
	    segs <- segmentations(object)
	    cat("number of segmentations calculated : ",nrow(segs),sep="")
	    cat('\n',sep="")	    
	    cat("least number of changepoints  = ",min(segs$m), " : maximum number of changepoints = ",max(segs$m),sep="")
	    cat('\n',sep="")
	    cat("number of iterations = ",object@iterations,sep="")
	    cat('\n',sep="")
            invisible()
})