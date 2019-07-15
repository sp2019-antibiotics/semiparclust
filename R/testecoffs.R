#' Validation Of ECOFF Using Only Scrucca's Model
#'
#' Tests, if ECOFF can be calculated for all observed data points using only Scrucca's model. Additionally it compares the ECOFF returned by fitting_fun to the official ECOFF values using a RSME.
#'
#' @param ZD data
#'
#' @return
#' \itemize{
#'  \item{errors_scrucca: number of ecoffs of observations with more than 4 components that could not be calculated}
#'  \item{RMSE_scrucca: root squared error}
#'}



testecoffs<- function(ZD){
    ecoff_only_scrucca = sapply(1:nrow(ZD), function(i) {print(i); tryCatch(fitting_fun(ZD[i, 4:48], fit_one_mode = F), error=function(error) NA)})
  (errors_scrucca<-nrow(ZD) - sum(is.na(ecoff_only_scrucca)))  # how many observations with more than 4 (i.e. 5 or more) components can not be estimated by scrucca


  ecoff_only_scrucca[ecoff_only_scrucca == 0] <- NA
  RMSE_scrucca = sqrt(mean((ecoff_only_scrucca[!(is.na(ZD$ECOFF) | is.na(ecoff_only_scrucca))] - ZD$ECOFF[!(is.na(ZD$ECOFF) | is.na(ecoff_only_scrucca))])^2))
  RMSE_scrucca
}
