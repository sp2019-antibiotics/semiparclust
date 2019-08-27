#' EM Density And ECOFF Value
#'
#' Plots the EM density and the ECOFF Value. Additionally it gives information on the number of components.
#'
#' @param ZD data
#' @param ant Antimicrobial
#' @param bac Bacterium
#'
#' @return plot of EM density and ECOFF value
#' @export
#'
#' @examples
#' data("ZD", package = "EUCASTData")
#' plotdens(ZD, "Ampicillin", "Escherichia coli")
#' plotdens(ZD, "Piperacillin", "Escherichia coli")
#' plotdens(ZD, "Mecillinam", "Escherichia coli")
#' @import graphics

plotdens<-function(ZD, ant, bac){
  binned_data1 <-  ZD[ZD$Antimicrobial == ant & ZD$Bacterium == bac, grepl("^Z", colnames(ZD)), drop = FALSE]
  set.seed(10) # because of the kmeans
  #final results:
  fitting_fun(binned_data1,fit_one_mode = T, plotting=T)
}
