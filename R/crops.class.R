setOldClass("memoised")
setOldClass("set")

#' An S4 class to encapsulate results produced by running a \code{crops} analysis
#'
#' @slot method A memoised version of the segmentation function provided by the user
#' @slot betas A set containing the beta values for which the sementation method was called during the \code{crops} analysis. These can be used in conjuction with \code{method} to efficiently recall the individual changepoint analysis results produced during the \code{crops} analysis 
#' @export
.crops.class<-setClass("crops.class",representation(method="memoised",betas="set"))
crops.class<-function(method,betas)
{
    .crops.class(method=method,betas=betas)
}





