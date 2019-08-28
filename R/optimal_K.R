#' optimal_K
#'
#' Assesses the optimal number of components GMM.
#'
#' @param data binned data
#' @param max_comp the maximal K which should be tried, function will try 1:max comp values of
#' @param criterion which criterion should be used, one can choose either AIC or BIC
#' @param take_resist logical, whether resistant observations - first column of the data - should be used, default is FALSE
#' @return
#' \itemize{
#'  \item{opt: optimal K - number of components which minimize chosen information criterion for the given data set}
#'  }
#' @export
#'
optimal_K <- function(data, max_comp, criterion, take_resist = F){
  crit_val <- sapply(1:max_comp, function(i) tryCatch(eval(parse(text=(paste("EM_fit(data, i , take_resist = take_resist, plot = F)$",criterion,sep="")))), error = function(error) NA))
  opt = which.min(crit_val)
  return(opt)
}
