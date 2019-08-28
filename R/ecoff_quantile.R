#' ECOFF Based On Quantile
#'
#' Finds ECOFF value based on the quantile of the GMM chosen by some information criterion.
#'
#' @param my_k  number of components of the GMM, can be calculated using optimal_K(data, round(sqrt(sum(binned_data!=0))), "BIC", take_resist = F)
#' @param res  result from the function EM_fit, can be calculated using EM_fit(binned_data, my_k, take_resist = F)
#' @param quant  quantile to use as a decision for the ECOFF value, default is 0.01
#' @return
#' \itemize{
#'  \item{out: ECOFF}
#' }
#' @export

ecoff_quantile = function(my_k, res, quant=0.01){
  vector <- sapply(seq(0,60,by=0.1), function(x) sum(sapply(1:my_k, function(k) res$pi[k] * pnorm(x, mean = res$mu[k], sd = res$sigma))))
  out <- seq(0,60,by=0.1)[which.min(abs(vector-quant))]
  return(out)
}
