#' Validation Of ECOFF
#'
#' Tests, if ECOFF can be calculated for all observed data points. Additionally it compares the ECOFF returned by fitting_fun to the official ECOFF values using a RSME.
#'
#' @param ZD data
#'
#' @return
#' \itemize{
#'  \item{errors: number of ecoffs that could not be calculated}
#'  \item{RMSE: root mean squared error}
#'  }
#' number of ecoffs that could not be calculated, RMSE

testecoff<-function(ZD){ # trycatch function catch the error
  new_ecoff = sapply(1:nrow(ZD), function(i) {print(i); tryCatch(fitting_fun(ZD[i, 4:48], fit_one_mode = T), error=function(error) NA)})   # print i for my debugging

  (errors<- sum(is.na(new_ecoff))) ## if this is equal to 0, we do not have any error.

  #excluding those ecoffs which are set to 0 (i.e. datasets with less observations than 5)
  new_ecoff[new_ecoff == 0] <- NA
  #root mean squared error for datasets which have some predefined ECOFF value in the dataset
  RMSE = sqrt(mean((new_ecoff[!(is.na(ZD$ECOFF) | is.na(new_ecoff))] - ZD$ECOFF[!(is.na(ZD$ECOFF) | is.na(new_ecoff))])^2))
  RMSE
}

