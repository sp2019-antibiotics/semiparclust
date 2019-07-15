#' Data Plots
#'
#' Plots given data.
#'
#' @param ZD data
#' @param anti Antimicrobial
#' @param bac Bacterium
#' @examples
#' plotdata(ZD, "Ampicillin", "Escherichia coli")
#' plotdata(ZD, "Piperacillin", "Escherichia coli")
#' plotdata(ZD, "Mecillinam", "Escherichia coli")
#' @import graphics
#' @export
plotdata <-function(ZD, anti, bac){
  ZDs <- ZD[ZD$Antimicrobial == anti & ZD$Bacterium == bac, grepl("^Z", colnames(ZD)), drop = FALSE]
  example <- data.frame(ZD = as.integer(gsub("^Z", "", colnames(ZDs))),
    Freq = unname(unlist(ZDs)))
  plot(Freq ~ ZD, data = example, type = "h")
}

#}
