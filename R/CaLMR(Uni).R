#' Univariable (single-exposure) Causal analysis of Latent exposures using Mendelian Randomization CaLMR (Uni)
#'
#' Causal analysis of Latent exposures using Mendelian Randomization (CaLMR) is an MR method that tests the causal relationships between
#' the outcome and the latent exposure using GWAS summary-level association statistics.
#' This function conducts CaLMR(Uni) test, assuming there is one latent exposure.
#' It is built under a two-sample MR framework and conducts Bayesian modeling using conjugate priors and Regression with Summary Statistics (RSS) Likelihood.
#'
#'@param sumtable a M*(K+1) data frame containing the GWAS summary data for K observable traits and the outcome.
#'@param Corr.mat a (K+1)*(K+1) estimated correlation matrix of the GWAS summary statistics.
#'                The 1st-Kth variables are related to the K observable traits, and the last variable corresponds to the outcome.
#'                *The order of the observable traits should match with the order in the 'traitvec' vector.
#'@param K number of observable traits
#'@param traitvec a vector containing the names of the observable traits. This should match with the column names in sumtable.
#'@param outcome the name of the outcome Y, and this should match with the column name in sumtable.
#'@param sign  a vector of length K, where the k-th element gives the pre-known sign of theta_k for the k-th trait in \code{traitvec}.
#'            Use 1 for positive, -1 for negative, and 0 if the sign is unknown.
#'@param T total number of iterations for the Gibbs sampler, with default T=3000
#'@param burnin length of burn-in period, with default burnin=1500.
#'@return
#'\itemize{
#' \item \code{calmr_result}: a vector contains:
#'           \code{calmr.p} (p-value for the causal effect theta),
#'           \code{calmr.rej} (testing for theta), and
#'           \code{calmr.sign} (sign of theta: 1, -1, or 0).
#' \item \code{CI}: a matrix of posterior 95% credible intervals for all model parameters after burn-in.
#' \item \code{mcmc.detail}: full MCMC chain for all parameters, can be used to check for convergence.
#'  }
#' @export


##############################################################
calmr_uni <- function(sumtable, Corr.mat, T=3000, burnin=1500, K = K, traitvec, outcome, sign) {

  if (K < 3) {
    stop("You must use at least 3 biomarkers.")
  }

  # Ensure Corr.mat has the same row and column orders as traitvec
  Corr.mat <- Corr.mat[traitvec, traitvec]

  if (nrow(Corr.mat) == K) {
    Corr.mat <- as.data.frame(Corr.mat)
    Corr.mat[[outcome]] <- 0
    Corr.mat <- rbind(Corr.mat, setNames(as.list(rep(0, ncol(Corr.mat))), colnames(Corr.mat)))
    rownames(Corr.mat)[nrow(Corr.mat)] <- outcome
    Corr.mat <- as.matrix(Corr.mat)
  }

  ########################################################################
  ## GWAS summary statistics
  sbetay  <- sumtable[, paste0("beta.", outcome)]
  ssigmay <- sumtable[, paste0("se.",   outcome)]
  sbetaBk <- ssigmaBk <- list()
  for (i in 1:K) {
    sbetaBk[[i]]  <- sumtable[, paste0("beta.", traitvec[i])]
    ssigmaBk[[i]] <- sumtable[, paste0("se.",   traitvec[i])]
  }

  M      <- length(sbetay)
  tau_X2 <- 6.666667e-05

  ########################################################################
  # Per-SNP summary-stat matrices in (M x (K+1)) layout
  S_mat <- matrix(0, nrow = M, ncol = K + 1)
  for (i in 1:K) S_mat[, i] <- sbetaBk[[i]]
  S_mat[, K + 1] <- sbetay

  se_mat <- matrix(0, nrow = M, ncol = K + 1)
  for (i in 1:K) se_mat[, i] <- ssigmaBk[[i]]
  se_mat[, K + 1] <- ssigmay

  # Per-SNP (K+1)x(K+1) covariance Omega_j (used in beta_X update),
  # per-SNP eta_sigma inverse (used in eta_theta update, cross-terms NOT scaled
  # by SEs to match the original Kronecker-product behavior),
  # and KxK biomarker-only subblock inverse (used in gamma update).
  Omega_list         <- vector("list", M)
  eta_sigma_inv_list <- vector("list", M)
  Omegak_inv_list    <- vector("list", M)

  for (j in 1:M) {
    Omega_j <- matrix(0, K + 1, K + 1)
    for (k in 1:(K + 1)) {
      for (i in 1:(K + 1)) {
        Omega_j[k, i] <- Corr.mat[k, i] * se_mat[j, k] * se_mat[j, i]
      }
    }
    Omega_list[[j]] <- Omega_j

    eta_sigma_j <- Corr.mat
    for (k in 1:K) {
      for (i in 1:K) {
        eta_sigma_j[k, i] <- Corr.mat[k, i] * se_mat[j, k] * se_mat[j, i]
      }
    }
    eta_sigma_j[K + 1, K + 1] <- ssigmay[j]^2
    eta_sigma_inv_list[[j]] <- solve(eta_sigma_j)

    Omegak_inv_list[[j]] <- solve(Omega_j[1:K, 1:K])
  }

  ########################################################################
  # Parameter initialization
  eta_theta <- c(0.5 * sign, 0.5)

  gamma_mat <- matrix(0.1, nrow = M, ncol = K)

  tau_k2 <- rep(1e-04, K)

  beta_X <- numeric(M)

  prior_alphak <- rep(0.0001, K)
  prior_betak  <- rep(0.0001, K)

  eta_simulated <- matrix(nrow = T, ncol = 2 * K + 1)

  ########################################################################
  # MCMC loop
  for (t in 1:T) {

    ############################
    ## A_theta is block-diagonal with eta per SNP, so t(A)%*%A = (eta^T eta) * I_M.
    ## C is diagonal -> each beta_X[j] is an independent rnorm draw.
    eta_vec <- eta_theta
    ata     <- sum(eta_vec^2)

    for (j in 1:M) {
      C_jj     <- as.numeric(crossprod(eta_vec, Omega_list[[j]] %*% eta_vec)) / ata^2
      resid_j  <- S_mat[j, ] - c(gamma_mat[j, ], 0)
      mu_raw_j <- sum(eta_vec * resid_j) / ata

      sig2_j  <- 1 / (1/C_jj + 1/tau_X2)
      mu_bx_j <- sig2_j * (1/C_jj) * mu_raw_j
      beta_X[j] <- rnorm(1, mean = mu_bx_j, sd = sqrt(sig2_j))
    }

    #############################
    tau_inv_diag <- diag(1/tau_k2, nrow = K)
    for (j in 1:M) {
      Omegak_inv_j <- Omegak_inv_list[[j]]
      Sigma_gk     <- solve(Omegak_inv_j + tau_inv_diag)
      resid_k      <- S_mat[j, 1:K] - eta_theta[1:K] * beta_X[j]
      mu_gk        <- Sigma_gk %*% (Omegak_inv_j %*% resid_k)
      gamma_mat[j, ] <- mvrnorm(n = 1, mu = mu_gk, Sigma = Sigma_gk)
    }

    #############################
    for (k in 1:K) {
      alpha_posterior_k <- prior_alphak[k] + M/2
      beta_posterior_k  <- prior_betak[k]  + 0.5 * sum(gamma_mat[, k]^2)
      tau_k2[k] <- 1 / rgamma(1, shape = alpha_posterior_k, rate = beta_posterior_k)
    }

    ##############################
    ## solve(eta_sigma / bx^2) = bx^2 * eta_sigma_inv  (avoids per-SNP inversion).
    sum_precision <- matrix(0, K + 1, K + 1)
    sum_weighted  <- numeric(K + 1)
    for (j in 1:M) {
      bx2        <- beta_X[j]^2
      scaled_inv <- bx2 * eta_sigma_inv_list[[j]]
      sum_precision <- sum_precision + scaled_inv

      S_j     <- S_mat[j, ]
      gamma_j <- c(gamma_mat[j, ], 0)
      sum_weighted <- sum_weighted +
        as.numeric(scaled_inv %*% ((S_j - gamma_j) / beta_X[j]))
    }
    Sigma_eta <- solve(sum_precision)
    mu_eta    <- Sigma_eta %*% sum_weighted
    eta_theta <- mvrnorm(n = 1, mu = mu_eta, Sigma = Sigma_eta)

    eta_simulated[t, ] <- c(eta_theta[K+1], eta_theta[1:K], tau_k2)
    if (t %% 10 == 0) cat("Iteration No: ", t, "\n")
  }

  ########################################################################
  eta_simulated <- as.data.frame(eta_simulated)
  ci <- t(sapply(1:ncol(eta_simulated), function(x){
    quantile(eta_simulated[-(1:burnin), x], c(0.025, 0.975))
  }))

  res.p <- rep(0, K + 1)
  for (i in 1:(K + 1)) {
    se      <- (ci[i, 2] - ci[i, 1]) / (2 * qnorm(0.975))
    z       <- mean(eta_simulated[-(1:burnin), i]) / se
    tmp.p   <- pnorm(-abs(z), log.p = TRUE)
    tmp.p   <- 2 * exp(tmp.p)
    res.p[i] <- tmp.p
  }
  if (sum(res.p[2:(K+1)] < 0.3) < 1) {
    stop("No significant thetak. Recommend to change biomarkers.")
  }

  idx <- which.min(res.p[2:(K+1)])
  if (mean(eta_simulated[-(1:burnin), (idx + 1)]) * sign[idx] < 0) {
    eta_simulated[, 1:(K+1)] <- eta_simulated[, 1:(K+1)] * (-1)
  }

  ci <- t(sapply(1:ncol(eta_simulated), function(x){
    quantile(eta_simulated[-(1:burnin), x], c(0.025, 0.975))
  }))
  bayes.rej <- 1
  if (ci[1, 1] < 0 & ci[1, 2] > 0) bayes.rej <- 0
  colnames(eta_simulated) <- rownames(ci) <- c("theta",
                                               paste0("theta", 1:K),
                                               paste0("tauk_", 1:K, "2"))

  bayes.sign <- 0
  if (mean(eta_simulated[-(1:burnin), 1]) > 0 & bayes.rej == 1) bayes.sign <-  1
  if (mean(eta_simulated[-(1:burnin), 1]) < 0 & bayes.rej == 1) bayes.sign <- -1

  Res <- list(calmr_result = c(calmr.p = res.p[1],
                               calmr.rej = bayes.rej,
                               calmr.sign = bayes.sign),
              CI = ci,
              mcmc.detail = eta_simulated)
  return(Res)
}
