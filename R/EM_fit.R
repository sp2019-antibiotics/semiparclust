#' EM Fit Of K-mixture Components
#'
#' For the estimation of GMM with only K=1 component we have used the function MASS::fitdistr().For the estimation of GMM with k=2,...,K components an EM algorithm is used to obtain needed output.
#'
#' @param data vector of binned data
#' @param K number of gaussian components for mixture models - at least two
#' @param take_resist logical, whether or not to consider first bin refering to resisant observations, default is F
#' @param plot logical, whether or not to plot final fit of GMM on original data, default is F
#' @param tolerance  tolerance for the Aitken-acceleration based stopping criterion, default is 0.001
#'
#' @return
#' \itemize{
#'  \item {pi: mixing parameters pik of final GMM}
#'  \item {mu: means muk of the final GMM}
#'  \item {sigma: standard deviation sigma of the final homoscedastic GMM}
#'  \item {loglik: vector of log-likelihood during the iterative procedure - in case of K=1 only last value of log-likelihood will be returned}
#'  \item {AIC: AIC of the GMM fitted to the given data}
#'  \item {BIC: BIC of the GMM fitted to the given data}
#'  \item {data_used: data which were used for the fitting procedure}
#'  \item {reconstr_data: reconstructed data of the given binned data}
#'  \item {bins: bins}
#' }
#' @export


EM_fit <- function(data, K, take_resist = F, plot = F, tolerance = 0.001){

  iter = 100

  if(take_resist == T){
    mydata <- as.numeric(data)
  }else{
    mydata <- as.numeric(data)[-1]
  }
  bins <- (50-length(mydata) + 1):50

  #reconstructed data
  reconstructed_data <- c()
  for (i in 1:length(mydata)){
    reconstructed_data = c(reconstructed_data, rep(bins[i], mydata[i]))
  }
  #total number of observations
  N <- length(reconstructed_data)

  if(K==1){
    my_dist = MASS::fitdistr(reconstructed_data, "normal")

    if(plot == T){
      hist(reconstructed_data, breaks = 5:max(bins), probability = T)
      lines( -100:100, dnorm(-100:100, my_dist$estimate[1], my_dist$estimate[2]))
    }

    return(list("pi" = c(1),"mu" =  my_dist$estimate[1],"sigma" = my_dist$estimate[2],
      "loglik" = my_dist$loglik, "AIC" =  AIC(my_dist), "BIC" = BIC(my_dist), "data_used" = mydata,"reconstr_data" = reconstructed_data,"bins" = bins))
  }else{


    #endpoints of bins
    b_j <- sapply(1:(length(mydata)-1), function(x) (bins[x] + bins[x + 1])/2)
    b_j[length(bins)] <- 60
    #starting points of bins
    a_j <- c(-10,b_j[-length(bins)])

    #starting values for EM obtained from kmeans clustering
    set.seed(666)
    kmeans_res <- suppressWarnings({kmeans(reconstructed_data, K, iter.max = 10, nstart = 10)})
    # we will not hear about not converging kmeans, we dont care anyway, as this is only for starting values
    mus <- sort(c(kmeans_res$centers))
    km_vars <- sapply(1:K, function(x) var(reconstructed_data[kmeans_res$cluster == x]))
    km_vars[is.na(km_vars)] <- 0 # if there is only 1 observation in a cluster, sgm would be inf: we set it to 0 for weighted mean
    weights <- c(sapply(1:K, function(x) sum(kmeans_res$cluster == x)))/N
    #starting value for sigma as weighted average over k-means cluster specific sigmas
    sgm = sqrt(weighted.mean(x = c(km_vars), w = weights))
    #equal starting probabilities for the mixture components
    pis <- rep(1/K, K)

    #empty matrices where values of parameters in each iteration will be recorded
    pi <- matrix(NA, K, iter+1)
    pi[,1] <- pis #first value is the initial from kmeans
    mu <- matrix(NA, K, iter+1)
    mu[,1] <- mus #first value is the initial from kmeans
    sigma <- c()
    sigma[1] <- sgm #first value is the initial from kmeans

    loglik <- c() #empty vector for storage of loglikelihood values over iterations

    #mixture density for K components
    dmix <- function(x){
      prob = matrix(sapply(1:K, function(k) (pi[k,i] * dnorm(x,mean=mu[k,i],sd=sigma[i]))), nrow = length(x))
      return (rowSums(prob))
    }

    difference = 100
    i = 1
    loglik_A = c()
    loglik_A[1] = -1000000*N #just a first value for stopping crit.

    #iterative procedure (formulae accordingly to report notes)
    while((difference>tolerance)&(i<iter)){
      m_j = mydata #binned data, in report as n_j

      H0_kj = sapply(1:K, function(k) {
        ifelse(b_j < mu[k, i],
          apply(pnorm(cbind(a_j, b_j),mean=mu[k,i],sd=sigma[i]), 1, diff),
          apply(pnorm(cbind(b_j, a_j),mean=mu[k,i],sd=sigma[i], lower.tail = FALSE), 1, diff))
      })
      p_j <- c(H0_kj%*%c(pi[,i]))
      H1_kj = sapply(1:K, function(k) dnorm(b_j,mean=mu[k,i],sd=sigma[i]) - dnorm(a_j,mean=mu[k,i],sd=sigma[i]))
      H2_kj = sapply(1:K, function(k) b_j*dnorm(b_j,mean=mu[k,i],sd=sigma[i]) - a_j*dnorm(a_j,mean=mu[k,i],sd=sigma[i]))

      E_tau = (H0_kj%*%diag(pi[,i])) / (p_j)
      E_tau_Y = ((H0_kj%*%diag(mu[,i]) - (sigma[i]^2)*H1_kj)%*%diag(pi[,i])) / (p_j)

      pi[,i+1]  <- (t(m_j)%*%E_tau)/sum(m_j) #new pi
      mu[,i+1] <- sort((t(m_j)%*%E_tau_Y)/(t(m_j)%*%E_tau)) #new mu
      E_tau_var = (((sigma[i]^2)*(H0_kj + H1_kj%*%(2*diag(mu[,i+1]) - diag(mu[,i])) - H2_kj) +H0_kj%*%diag(c(mu[,i+1] - mu[,i]))^2 )%*%diag(pi[,i]))/(p_j)

      sigmas <-  (t(m_j)%*%E_tau_var)/(t(m_j)%*%E_tau )
      sigma[i+1] <- sqrt((t(m_j)%*%E_tau %*% t(sigmas))/N)

      loglik[i] <- sum(m_j*log(p_j))

      if(i>2){ # computing of aitken acceleration based stopping criterion
        a = (loglik[i] - loglik[i-1])/(loglik[i-1] - loglik[i-2])
        loglik_A[i-1] = loglik_A[i-2] + (1/(1-a))*(loglik[i] - loglik[i-1])
        difference = abs(loglik_A[i-1]-loglik_A[i-2])
      }

      i=i+1
    }

    #together we have 2*K parameters (K-1 for pi, K for mu, 1 for sigma)
    AIC = 2*2*K - 2*loglik[i-1]
    BIC = log(N)*2*K - 2*loglik[i-1]

    final_density <- function(x){
      prob = matrix(sapply(1:K, function(k) (pi[k,i] * dnorm(x,mean=mu[k,i],sd=sigma[i]))), nrow = length(x))
      return (rowSums(prob))
    }

    #plot it if needed
    if(plot == T){
      hist(reconstructed_data, breaks = 5:max(bins), probability = T, col = "red",
        main = "EM density and ECOFF value",
        xlab = "ZD")
      lines( seq(0,50,0.001), final_density(seq(0,50,0.001)), lwd = 2)
    }


    return(list("pi" = pi[,i],"mu" =  mu[,i],"sigma" = sigma[length(sigma)],
      "loglik" = loglik, "AIC" =  AIC, "BIC" = BIC, "data_used" = mydata,"reconstr_data" = reconstructed_data,"bins" = bins))
  }
}

