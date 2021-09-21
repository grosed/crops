setOldClass("set")
.crops.class<-setClass("crops.class",representation(log="set"))
crops.class<-function(log)
{
    .crops.class(log=log)
}


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


setMethod("plot",signature=list("crops.class","data.frame"),
          function(x,y)
          {
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

setMethod("plot",signature=list("data.frame","crops.class"),
          function(x,y)
          {
	    return(plot(y,x))
          })


setMethod("plot",signature=list("crops.class","missing"),
          function(x)
          {
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