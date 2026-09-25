setOldClass("memoised")
setOldClass("set")
.crops.class<-setClass("crops.class",representation(method="memoised",betas="set"))
crops.class<-function(method,betas)
{
    .crops.class(method=method,betas=betas)
}





