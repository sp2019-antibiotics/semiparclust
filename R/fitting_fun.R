
#' Final ECOFF Value
#'
#' As the ECOFF value can not always be calculated through semi-parametric clustering, this function decides how the ECOFF should be calculated (semi_param_clust or ecoff_quantile) and returns the final ECOFF. If there is not enough information in the data, this function will return ECOFF to be 0. This can happen in the case when we have less than 5 observations(i.e. 4 and less).
#'
#' @param data  data which is used for analysis (one row of the data set)
#' @param fit_one_mode  logical, whether function ecoff quantile should be called for the data, if they can not be clustered to 2 clusters by semi_param_clust, default is F
#' @param quant quantile which should be used for the function ecoff_quantile, default is 0.01
#' @param take_resist logical, whether the resistant - in our case Z6, i.e. the first column - should be used for the fitting and further procedure, default is F
#' @param plotting  logical, whether the histogram of the data with fitted GMM and ECOFF value should be plotted, default is F
#' @return
#' \itemize{
#'  \item{my_ecoff: ECOFF}
#' }
#' @export
#' @import stats


fitting_fun <- function(data, fit_one_mode=NULL, quant=0.01, take_resist=F, plotting = F){
  #data <- data[1, grepl("^Z", colnames(data))]
  # data: one row only
  binned_data <- as.numeric(data)

  n_nonzero = ifelse(take_resist==F, sum(binned_data[-1]!=0), sum(binned_data!=0))
  n_obs = ifelse(take_resist==F, sum(binned_data[-1]), sum(binned_data))
  if((n_obs<5)&(n_nonzero<2)){
    # if the number of observations is lower than 5, or number of non-zero bins is lower than 2, we set ecoff value to 0
    return("my_ecoff" = 0)
  }

  max_components = ifelse(take_resist==F, round(sum(binned_data[-1]!=0)/2), round(sum(binned_data)/2))

  #optimal number of components
  my_k <- optimal_K(binned_data, max_components, "BIC", take_resist = take_resist) ## must be BIC !!
  #result fit
  res <- EM_fit(binned_data, my_k, take_resist = take_resist)
  semi_res = semi_param_clust(res, my_k)$ecoff

  if((fit_one_mode==T)&(is.na(semi_res))){
    semi_res = ecoff_quantile(my_k=my_k, res=res, quant=quant)
  }


  my_dens <- function(x){
    prob = matrix(sapply(1:my_k, function(k) (res$pi[k] * dnorm(x,mean=res$mu[k],sd=res$sigma))), nrow = length(x))
    return (rowSums(prob))
  }

  if(plotting == T){
    hist(as.numeric(res$reconstr_data), breaks = 5:max(res$bins),
      probability = T, col = "red", main = "EM density and ECOFF value",
      xlab = "ZD")
    lines( seq(0,50,0.01), my_dens(seq(0,50,0.01)), lwd = 2)
    abline(v = semi_res, col = "darkblue", lty = 2, lwd = 2)
    legend("topright",
      legend = c(paste("Components:", as.character(my_k)),
        paste("ECOFF:", as.character(semi_res))),
      bty = "n")
  }

  return("my_ecoff" = semi_res)

}
