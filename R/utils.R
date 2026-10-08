
check_crops_arguments <- function(method,beta_min,beta_max)
{
    if(!is(method,'function'))
    {
      stop("method should be a function. See crops documentation")
    }
    if(!is(beta_min,'numeric') | length(beta_min) > 1 | beta_min < 0.0)
    {
      stop("beta_min should be a single positive numeric value. See crops documentation")
    }
    if(!is(beta_max,'numeric') | length(beta_min) > 1 | beta_min < 0.0)
    {
      stop("beta_max should be a single positive numeric value. See crops documentation")
    }
    if(beta_min >= beta_max)
    {
       stop("beta_max should be greater than beta_min. See crops documentation.")
    }
    return()
}


check_method_return_values <- function(res)
{
    if(!is(res,'list'))
     {
       stop("method should be a function returning a list. See crops documentation.")
     }
    if(!is(res[[1]],'numeric') | length(res[[1]]) > 1 | res[[1]] < 0.0) 
    {
      msg <- paste("First value in the list returned by method should be a single positive real value. See crops documentation.")
      msg <- paste(msg,'\n')
      stop(msg)
    }
    if(!is(res[[2]],'numeric'))
    {
      msg <- paste("Second value in the list returned by method should be a single integer numeric value indicating the number of changepoints. See crops documentation.")
      stop(msg)
    }
    return()
}