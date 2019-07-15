#' Semi-Parametric Clustering According To Scrucca
#'
#' Uses Semi-Parametric Clustering to part the data into two clusters. The border between those two clusters is estimated to be the ECOFF.
#'
#' @param res  result from the function EM fit
#' @param my_k number of components of the GMM
#' @param plot if TRUE, the function will return a result nd plot the fitted density and the mode of the function - set FALSE by default
#'
#' @return
#' \itemize{
#'  \item{Cluster1}
#'  \item{Cluster2}
#'  \item{ecoff}
#'  }

#' @references Scrucca (2016) <doi:10.1016/j.csda.2015.01.006>
semi_param_clust <- function(res, my_k, plot = F){
  # res: result from EM
  # my_k: optimal K


  #final density of GMM
  my_dens <- function(x){
    prob = matrix(sapply(1:my_k, function(k) (res$pi[k] * dnorm(x,mean=res$mu[k],sd=res$sigma))), nrow = length(x))
    return (rowSums(prob))
  }

  unbinned_data <- res$reconstr_data
  N = length(unbinned_data)
  binned_data <- res$data_used
  bins <- res$bins

  #fitted density
  #plot it if needed
  if(plot == T){
    plot(seq(0,50,length.out = 10000), my_dens(seq(0,50,length.out = 10000)), type = "l",
      xlab = "x", ylab = "density(x)")
  }


  #number of grid points for p (used for estimating sample level set)
  n_grid <- 50

  x <- seq(0,60,length.out = 10000)
  fx <- my_dens(x)
  myc <- seq(max(fx), 0, length.out = n_grid)

  m <- sapply(1:length(myc), function(p) sum(rle(fx > myc[p])$values == T))


  #Delauney triangulation (for binned data with integer-bins) is obtained by rounding the values to the closest integer
  #Connected sets are calculated for each p (find observations from High Density Interval corresponding to the value p)

  #mode function for different p
  #plot it if needed
  if(plot == T){
    plot(myc, m, xlim = rev(range(myc)),xlab =  "c = f(x)" , ylab = "m(c)", main = "Mode Function","s", lwd = 2)
  }


  #number of cores is equal to the number of increments of the mode function
  n_cores <- sum(diff(m) > 0)

  #For our project we consider always 2 clusters to identify ECOFF value.
  #For the value of p we choose the mean value of p corresponding to the longest sequence of 2s in the mode function
  #If there are more sequences of 2 of the same length, then we choose the first one
  first2 = sum(rle(m)$lengths[1:(which.max((rle(m)$values==2)*(rle(m)$lengths))-1)])+1 #index of first p where m(p) = 2 in the suitable sequence
  last2 = sum(rle(m)$lengths[1:(which.max((rle(m)$values==2)*(rle(m)$lengths)))]) #index of last p where m(p) = 2 in the suitable sequence

  start_c <- mean(myc[first2:last2]) #p for the first step of cluster-assignment

  if(sum(rle(fx > start_c)$values == T) == 1){
    return(list("ecoff" = NA))
  }else{

    core1 <- round(x[rle(fx > start_c)$lengths[1] + 1]) : round(x[sum(rle(fx > start_c)$lengths[1:2])]) #round the values to obtain the bins in the connected set of the first cluster
    core2 <-  round(x[sum(rle(fx > start_c)$lengths[1:3]) + 1]) : round(x[sum(rle(fx > start_c)$lengths[1:4])]) #round the values to obtain the bins in the connected set of the second cluster

    if(max(core1)==min(core2)){ return(list("ecoff" = NA))}

    #alocated binned data
    binned_aloc_data <- binned_data[which(bins %in% c(core1, core2))]

    #alocated bins
    bins_aloc <- c(core1, core2)

    #unalocated data
    unbinned_unaloc <- unbinned_data[!unbinned_data %in% bins_aloc]

    #cores of clusters
    data_core1 <- unbinned_data[unbinned_data %in% core1]
    data_core2 <- unbinned_data[unbinned_data %in% core2]

    if(length(data_core1)*length(data_core2)==0){
      return(list("ecoff" = NA))
    }

    N_aloc <- length(c(data_core1, data_core2)) #number of alocated observations
    N_unaloc <- N - N_aloc
    while(N_unaloc > 0){
      ##fit GMM on alocated data
      pi_aloc <- c(length(data_core1), length(data_core2))/N_aloc #mixing prob
      mu_aloc <- c(mean(data_core1), mean(data_core2)) #mu values
      sd_datacore1 = ifelse(is.na(sd(data_core1)), 0, sd(data_core1))
      sd_datacore2 = ifelse(is.na(sd(data_core2)), 0, sd(data_core2))
      sd_aloc <- weighted.mean(c(sd_datacore1, sd_datacore2), w = pi_aloc)  #sigma is the same for both clusters
      #calculating log/ratio (risk of each unalocated observation) of being in certain cluster

      log_ratio <- function(x){
        z1 = (pi_aloc[1]*dnorm(x, mean = mu_aloc[1], sd =  sd_aloc))/
          (1 - pi_aloc[1]*dnorm(x, mean = mu_aloc[1], sd =  sd_aloc))
        z2 = (pi_aloc[2]*dnorm(x, mean = mu_aloc[2], sd = sd_aloc))/
          (1 - pi_aloc[2]*dnorm(x, mean = mu_aloc[2], sd = sd_aloc))
        r1 = log(z1)
        r2 = log(z2)
        return(list("r1" = r1, "r2" = r2))
      }

      #log-ratios
      lr_res <- log_ratio(unbinned_unaloc)

      #threshold value (of empirical quantile) for allocation
      q = sqrt(N_aloc/N)

      #empirical quantiles of log-ratios for each cluster
      quantile(lr_res$r1, probs = q)
      quantile(lr_res$r2, probs = q)

      #candidates for alocation
      ind1 <- which(lr_res$r1 >= quantile(lr_res$r1, probs = q))
      ind2 <- which(lr_res$r2 >= quantile(lr_res$r2, probs = q))

      #if one observation has log-ratio higher than threshold for both clusters, then we choose
      #the cluster having higher log-ratio (i.e. higher chance of observation being in this cluster)
      to1 = unique(unbinned_unaloc[ind1])
      to2 = unique(unbinned_unaloc[ind2])
      is.higher = sapply(1:length(unbinned_unaloc), function(r) which.max(c(lr_res$r1[r],lr_res$r2[r])))
      to.which = ((lr_res$r1 >= quantile(lr_res$r1, probs = q))==(lr_res$r2 >= quantile(lr_res$r2, probs = q)))*(lr_res$r1 >= quantile(lr_res$r1, probs = q))*is.higher
      ind1 = ind1[!ind1 %in% which(to.which==2)]
      ind2 = ind2[!ind2 %in% which(to.which==1)]

      add_to1 <- unbinned_unaloc[ind1]
      add_to2 <- unbinned_unaloc[ind2]

      # adding new datapoints to clusters
      data_core1 <- sort(c(data_core1, add_to1))
      data_core2 <- sort(c(data_core2, add_to2))

      # removing newly allocated from unallocated
      unbinned_unaloc <- unbinned_unaloc[!unbinned_unaloc %in% c(add_to1,add_to2)]

      # adding new bins to cores
      core1 <- c(core1, unique(add_to1))
      core2 <- c(core2, unique(add_to2))

      N_unaloc = length(unbinned_unaloc)
      #alocated binned data
      binned_aloc_data <- binned_data[which(bins %in% c(core1, core2))]
      #alocated bins
      bins_aloc <- c(core1, core2)
      #unalocated data
      unbinned_unaloc <- unbinned_data[!unbinned_data %in% bins_aloc]
      #cores of clusters
      data_core1 <- unbinned_data[unbinned_data %in% core1]
      data_core2 <- unbinned_data[unbinned_data %in% core2]

      N_aloc <- length(c(data_core1, data_core2)) #number of alocated observations
    }

    #data from cluster 1
    data_core1
    #data from cluster 2
    data_core2

    #our ECOFF is the border value between last element of the first cluster and the first element of the second cluster
    ecoff = mean(c(max(data_core1),min(data_core2)))

    return(list("Cluster1" = data_core1, "Cluster2" = data_core2, "ecoff" = ecoff))
  }
}  # end of semi-parametric clustering
