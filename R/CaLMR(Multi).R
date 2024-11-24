#'Multivariable (multi-exposure) Causal analysis of Latent exposures using Mendelian Randomization CaLMR(Multi)
#'
#' Causal analysis of Latent exposures using Mendelian Randomization (CaLMR) is an MR method that tests the causal relationships between
#' the outcome and the latent exposures using GWAS summary-level association statistics.
#' This function conducts CaLMR(Multi) test, assuming there are multiple latent exposures.
#' It is built under a two-sample MR framework and conducts Bayesian modeling using conjugate priors and Regression with Summary Statistics (RSS) Likelihood.
#'
#'@param sumtable a M*(K+1) data frame containing the GWAS summary data for K observable traits and the outcome.
#'@param Corr.mat a (K+1)*(K+1) estimated correlation matrix of the GWAS summary statistics.
#'                The 1st-Kth variables are related to the K observable traits, and the last variable corresponds to the outcome.
#'                *The order of the observable traits should match with the order in the 'traitvec' vector.
#'@param grp a list of length L with the sub-list grp[[l]] containing the names of the observable traits associated with l-th latent exposure.
#'@param L number of latent exposures
#'@param K number of observable traits
#'@param traitvec a vector containing the names of the observable traits. This should match with the column names in sumtable.
#'@param outcome the name of the outcome Y, and this should match with the column name in sumtable.
#'@param sign a list of length L. Each sub-list sign[[l]] should be a vector of length K containing the pre-known signs of theta_kl.
#'             [1 means positive, -1 means negative, and 0 means this trait is not associated with l-th latent exposure or lack of pre-known signs]
#'@param T total number of iterations for the Gibbs sampler, with default T=3000
#'@param burnin length of burn-in period, with default burnin=1500.
#'@return a list of the testing results based on CaLMR
#'        - bayes.rej: a vector of length L showing the existences of causal effects for each latent exposure
#'        - bayes.sign: a vector of length L showing the the directions of the significant causal effects. 0 means no causal effects.
#'        - cor.X: estimated correlation matrix for the latent exposures
#'        - par: posterior samples
#'        - ci: 95% credible intervals for model parameters after dropping the burnin period.

##############################################################
calmr_multi <- function(sumtable, Corr.mat, grp, L, K, traitvec, outcome, sign, T=3000, burnin=1500) {
  #########################################################################
  orderref <- cbind(traitvec, num=1:K)
  grp.n <- grp
  for (l in 1:L){
    grp.n[[l]] <- sort(as.numeric(orderref[which(orderref[,1] %in% grp[[l]]), 2]))
  }
  grp <- grp.n

  ########################################################################
  ## Define the parameters coming from the summary data
  ##### GWAS summary statistics for Y (outcome):
  sbetay = sumtable[,paste0("beta.",outcome)]
  ssigmay = sumtable[,paste0("se.",outcome)]
  ##### GWAS summary statistics for biomarkers:
  sbetaBk = ssigmaBk = list()
  for (i in 1:K){
    sbetaBk[[i]] = sumtable[,paste0("beta.",traitvec[i])]
    ssigmaBk[[i]] = sumtable[,paste0("se.",traitvec[i])]
  }

  M <- length(sbetay)  # total number of IVs


  ########################################################################
  # The first K*M row is coef. for biomarkers, and the last M row is coef. for the outcome Y
  ## Matrix S is the summary coefficients, dimension = (K+1)*M by 1
  S <- matrix(0, nrow = (K+1)*M, ncol = 1)
  for (i in 1:K){
    for (j in 1:M){
      S[(i-1)*M+j,1]=sbetaBk[[i]][j]
    }
  }
  ## summary coef for the outcome Y
  S[(K*M+1):((K+1)*M),1]=sbetay

  ## The diagonal entries of the covariance matrix comes from std error in the summary data
  cov.mat <- Corr.mat %x% diag(M)
  ### Enter the variances of biomarkers
  for(i in 1:K){
    for (m in 1:K) {
      for (j in 1:M){
        cov.mat[(i-1)*M+j,(m-1)*M+j] = Corr.mat[i,m]*ssigmaBk[[i]][j]*ssigmaBk[[m]][j]
      }
    }
  }
  ### Enter the variances of the outcome Y
  for (j in 1:M){
    cov.mat[K*M+j,K*M+j]=ssigmay[j]^2
  }


  ########################################################################
  # Initiation of parameters (to be updated) in the model
  ## Vector of theta's: the first K elements are theta_k's
  ##                    the last element is theta, which is our main interest
  eta_theta = matrix(0, nrow=K+1, ncol=L) # initialize the theta's
  for (l in 1:L){
    eta_theta[grp[[l]],l] = 0.5  # distinguish the non-zero theta's
    eta_theta[K+1,l] = 0.5
  }

  ## B_gamma: dimension = (K+1)*M by 1
  B_gamma <- matrix(1e-2, nrow = (K+1)*M, ncol = 1)
  B_gamma[(K*M+1):((K+1)*M),]=0  ## FIXED: The last M rows are 0 (WILL NOT UPDATE)
  # initialize tau_k^2
  tau_k2 <- rep(1e-04,K)
  tau_X2 <- rep(6.666667e-05,L)
  D <- diag(sqrt(tau_X2))
  R <- matrix(0.5, ncol=L, nrow=L)
  diag(R) <- rep(1, L)
  H <- D%*%R%*%D
  beta_X <- matrix(1e-2, nrow=M, ncol=L)
  ########################################################################
  # Initiation of fixed parameters in the model
  # Define hyperparameters
  prior_alphak <- rep(0.0001,K); prior_betak <- rep(0.0001,K)  ##

  ########################################################################
  ########################################################################
  ########################################################################
  # Set up
  n.grp <- sum(sapply(grp, length))
  eta_simulated=matrix(nrow=T,ncol=n.grp+K+L)
  corr.X = matrix(nrow=T,ncol=L*(L-1)/2)

  # Start of iteration
  for (t in 1:T) {
    ############################
    # Update beta_X

    A_tem = solve(t(eta_theta)%*%eta_theta)%*%t(eta_theta)
    # update beta_x by SNPs
    for(j in 1:M){
      S_j=S[c(0:K)*M+j,]
      gamma_j=B_gamma[c(0:K)*M+j,]
      covmatj = matrix(0, K+1 ,K+1)
      for(i in 1:K){
        for(k in 1:K) {
          covmatj[i,k] = Corr.mat[i,k]*ssigmaBk[[i]][j]*ssigmaBk[[k]][j]
        }
      }
      covmatj[K+1, K+1] = ssigmay[j]^2
      Cj <- A_tem%*% covmatj %*% t(A_tem)
      sigma_betax <- solve(solve(Cj) + solve(H))
      mu_betax <- sigma_betax %*% (solve(Cj)%*%A_tem %*%(S_j-gamma_j) )

      beta_X[j,] <- mvrnorm(n =1, mu = mu_betax, Sigma = sigma_betax)
    }

    # update the between-exposure correlation matrix
    # this will be used to update the covariance matrix between exposures
    # but the diagonal entries of the covariance matrix are fixed at tau_X2
    W <- rinvwishart( nu=M+L+1, S=t(beta_X)%*%(beta_X)+diag(tau_X2))
    # H<-riwish(v=M+L+1,S=t(beta_X)%*%(beta_X)+diag(tau_X2))
    R <- solve(diag(sqrt(diag(W)))) %*% W %*% solve(diag(sqrt(diag(W))))
    H <- D%*%R%*%D
    corr.X[t,] <- R[upper.tri(R)]

    #############################
    # Update the first K*M entries of B_gamma
    # To save time, update gamma_k by SNPs

    for (j in 1:M) {
      Omegak = matrix(0, K, K)
      for (k in 1:K) {
        for (i in 1:K){
          Omegak[k, i] = Corr.mat[k,i]*ssigmaBk[[k]][j]*ssigmaBk[[i]][j]
        }
      }
      Omegak.solve = solve(Omegak)
      nus <- matrix(0, nrow=K, ncol=1)
      for (l in 1:L){
        nus <- nus + beta_X[j,l]*eta_theta[1:K,l]
      }
      Sigma_B_gamma_k <- solve(Omegak.solve+diag(1/tau_k2))
      mu_B_gamma_k <- Sigma_B_gamma_k %*% Omegak.solve%*%(S[c((0:(K-1))*M + j),]-nus)
      B_gamma[c((0:(K-1))*M + j),]<-mvrnorm(n=1, mu=mu_B_gamma_k, Sigma=Sigma_B_gamma_k)
    }


    #############################
    # Update tau_k2
    for (k in 1:K) {
      alpha_posterior_k = prior_alphak[k] + M/2
      beta_posterior_k = prior_betak[k] + 1/2*sum(B_gamma[((k-1)*M+1):(k*M),]^2)
      tau_k2[k] <- 1 / rgamma(1, shape = alpha_posterior_k, rate  = beta_posterior_k)
    }
    ##############################
    # Update eta_theta
    for (l in 1:L){
      index.l <- c(grp[[l]], K+1)
      kforl <- length(grp[[l]])
      theta_otherl <- c(1:L)[c(1:L)!=l]
      sum.betaX_etasigma=matrix(0,kforl+1,kforl+1)
      tmp.mu=matrix(0,kforl+1,1)

      for (j in 1:M){
        cov.mat_l = cov.mat[((index.l-1)*M+j), ((index.l-1)*M+j)]
        tmp.solve = solve(cov.mat_l/beta_X[j,l]^2)
        sum.betaX_etasigma = sum.betaX_etasigma + tmp.solve
        nus <- matrix(0, nrow=kforl+1, ncol=1)
        for (v in theta_otherl){
          nus <- nus + beta_X[j,v]*eta_theta[index.l,v]
        }
        S_j=S[c(index.l-1)*M+j,]
        gamma_j=B_gamma[c(index.l-1)*M+j,]
        tmp.mu = tmp.mu + tmp.solve %*% ((S_j-nus-gamma_j)/beta_X[j,l])
      }
      Sigma_eta <- solve(sum.betaX_etasigma)
      mu_eta <- Sigma_eta %*% tmp.mu
      eta_theta[index.l,l] <- mvrnorm(n = 1, mu = mu_eta, Sigma = Sigma_eta)
    }


    # save the results of each iteration
    eta_simulated[t,] <- c(eta_theta[K+1,], unlist(lapply(1:L, function(i) eta_theta[grp[[i]], i])), tau_k2)
    if (t %% 10 == 0) print(t) ##
  }

  colnames(eta_simulated) = c(paste0('theta', 1:L),
                              unlist(lapply(1:L, function(i) paste0('theta', grp[[i]], i))),
                              paste0('h2gamma', 1:K))

  eta_simulated<-as.data.frame(eta_simulated)

  for(l in 1:L){
    if(mean(eta_simulated[-(1:burnin),paste0("theta",grp[[l]][1], l)])*sign[[l]][grp[[l]][1]] < 0) {
      col<-c(paste0("theta",l), paste0('theta',grp[[l]],l))
      eta_simulated[,col]<-eta_simulated[,col]*(-1)
    }
  }

  # calculate the 95% CI
  ci <- t(sapply(1:ncol(eta_simulated), function(x){quantile(eta_simulated[-(1:burnin),x], c(0.025, 0.975))}))
  rownames(ci) <- colnames(eta_simulated)

  # get the estimated correlation matrix for the exposures
  corr_sub <- corr.X[-(1:burnin), ]
  if (is.vector(corr_sub)) {
    corr_sub <- matrix(corr_sub, ncol = 1)
  }
  cor.X <- diag(L)
  cor.X[upper.tri(cor.X)] <- colMeans(corr_sub)
  cor.X[lower.tri(cor.X)] <- t(cor.X)[lower.tri(cor.X)]

  # test result
  bayes.rej=rep(1, L)
  for(l in 1:L) {
    if(ci[paste0("theta",l),1]<0 & ci[paste0("theta",l),2]>0) {bayes.rej[l]=0}
  }
  # sign of thetas
  bayes.sign=rep(0, L)
  for(l in 1:L) {
    if(mean(eta_simulated[-(1:burnin),paste0("theta",l)])>0 & bayes.rej[l]==1) {bayes.sign[l]=1}
    if(mean(eta_simulated[-(1:burnin),paste0("theta",l)])<0 & bayes.rej[l]==1) {bayes.sign[l]=-1}
  }

  ## output
  Res <- list(bayes.rej = bayes.rej, bayes.sign = bayes.sign, cor.X = cor.X, par = eta_simulated, ci=ci)
  return(Res)
}

