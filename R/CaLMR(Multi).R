#' Multivariable (multiple exposure) Causal analysis of Latent exposures using Mendelian Randomization CaLMR (Multi)
#'
#' Causal analysis of Latent exposures using Mendelian Randomization (CaLMR) is an MR method that tests the causal relationships between
#' the outcome and multiple latent exposures using GWAS summary-level association statistics.
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
#'@param T total number of iterations for the Gibbs sampler, with default T=3000
#'@param burnin length of burn-in period, with default burnin=1500.
#'@param sign a list of length L. Each sub-list sign[[l]] should be a named or positional vector of length K,
#'            where the k-th element gives the pre-known sign of theta_kl for the k-th trait in \code{traitvec}.
#'            Use 1 for positive, -1 for negative, and 0 if the trait is not associated with the l-th latent exposure or if the sign is unknown.
#'            IMPORTANTLY, indexing follows the global position in \code{traitvec} (1 to K), not the local position within \code{grp[[l]]}.
#'@return
#'\itemize{
#' \item \code{calmr_result}: a data frame with one row per latent exposure with at least one significant corresponding theta_kl,containing three columns:
#'           \code{calmr.p} (p-value for the causal effect theta_l),
#'           \code{calmr.rej} (testing for this theta_l), and
#'           \code{calmr.sign} (sign of theta_l: 1, -1, or 0).
#'            Latent exposures with no significant theta_kl are excluded from this data frame.
#' \item \code{Res_latent_detail}: a list of length L. For latent exposure with at least one significant corresponding theta_kl,
#'          each entry is a list containing \code{calmr.p}, \code{calmr.rej}, \code{calmr.sign},
#'          and \code{thetakl.p} (a vector of p-values for each theta_kl).
#'          For latent exposures with no significant theta_kl, the entry is a character string recommending changing biomarkers.
#' \item \code{CI}: a matrix of posterior 95% credible intervals for all model parameters after burn-in.
#' \item \code{cor.X}: an L*L estimated correlation matrix among the L latent exposures.
#' \item \code{mcmc.detail}: full MCMC chain for all parameters, can be used to check for convergence.
#'}
#' @export

##############################################################
calmr_multi <- function(sumtable, Corr.mat, grp, L, K, T=3000, burnin=1500, traitvec, outcome, sign) {

  orderref <- cbind(traitvec, num = 1:K)
  grp.n <- grp
  for (l in 1:L) {
    grp.n[[l]] <- sort(as.numeric(orderref[which(orderref[, 1] %in% grp[[l]]), 2]))
  }
  grp <- grp.n

  ########################################################################
  ## GWAS summary statistics
  sbetay  <- sumtable[, paste0("beta.", outcome)]
  ssigmay <- sumtable[, paste0("se.",   outcome)]
  sbetaBk <- ssigmaBk <- list()
  for (i in 1:K) {
    sbetaBk[[i]]  <- sumtable[, paste0("beta.", traitvec[i])]
    ssigmaBk[[i]] <- sumtable[, paste0("se.",   traitvec[i])]
  }
  M <- length(sbetay)

  ########################################################################
  ## [Change 1] Pre-compute per-SNP summary stats and covariance inverses ONCE.
  ## Eliminates the (K+1)M x (K+1)M Kronecker product and ~M + L*M per-iter inversions.

  # Per-SNP summary-stat matrix in (M x (K+1)) layout
  S_mat <- matrix(0, nrow = M, ncol = K + 1)
  for (i in 1:K) S_mat[, i] <- sbetaBk[[i]]
  S_mat[, K + 1] <- sbetay

  # Per-SNP (K+1)x(K+1) cov for beta_X update (cross-terms = 0) and KxK biomarker
  # subblock inverse for gamma update.
  message("Pre-computing per-SNP covariance matrices...")
  betax_covmat_list <- vector("list", M)
  Omegak_inv_list   <- vector("list", M)

  for (j in 1:M) {
    covmatj <- matrix(0, K + 1, K + 1)
    for (i in 1:K) {
      for (k in 1:K) {
        covmatj[i, k] <- Corr.mat[i, k] * ssigmaBk[[i]][j] * ssigmaBk[[k]][j]
      }
    }
    covmatj[K + 1, K + 1] <- ssigmay[j]^2
    betax_covmat_list[[j]] <- covmatj
    Omegak_inv_list[[j]]   <- solve(covmatj[1:K, 1:K])
  }

  # Per-group, per-SNP sub-covariance inverses for eta_theta update.
  # eta_covmat differs from betax_covmat only in the biomarker-outcome cross-terms,
  # which stay as raw Corr.mat values (matching the original Kronecker behavior).
  index_lists <- lapply(1:L, function(l) c(grp[[l]], K + 1))

  message("Pre-computing per-group sub-covariance inverses...")
  eta_sub_inv <- vector("list", L)
  for (l in 1:L) {
    idx <- index_lists[[l]]
    eta_sub_inv[[l]] <- vector("list", M)
    for (j in 1:M) {
      eta_covmatj <- matrix(0, K + 1, K + 1)
      for (i in 1:K) {
        for (k in 1:K) {
          eta_covmatj[i, k] <- Corr.mat[i, k] * ssigmaBk[[i]][j] * ssigmaBk[[k]][j]
        }
      }
      eta_covmatj[K + 1, K + 1] <- ssigmay[j]^2
      for (i in 1:K) {
        eta_covmatj[i, K + 1] <- Corr.mat[i, K + 1]
        eta_covmatj[K + 1, i] <- Corr.mat[K + 1, i]
      }
      eta_sub_inv[[l]][[j]] <- solve(eta_covmatj[idx, idx])
    }
  }
  message("Pre-computation complete.")

  ########################################################################
  # Parameter initialization
  eta_theta <- matrix(0, nrow = K + 1, ncol = L)
  for (l in 1:L) {
    eta_theta[grp[[l]], l] <- 0.5
    eta_theta[K + 1, l]    <- 0.5
  }

  ## [Change 2] gamma stored as M x K matrix (drops the always-zero outcome rows)
  gamma_mat <- matrix(1e-2, nrow = M, ncol = K)

  tau_k2 <- rep(1e-04, K)
  tau_X2 <- rep(6.666667e-05, L)
  D <- diag(sqrt(tau_X2))
  R <- matrix(0.5, ncol = L, nrow = L); diag(R) <- 1
  H <- D %*% R %*% D
  beta_X <- matrix(1e-2, nrow = M, ncol = L)

  prior_alphak <- rep(0.0001, K)
  prior_betak  <- rep(0.0001, K)

  n.grp <- sum(sapply(grp, length))
  eta_simulated <- matrix(nrow = T, ncol = n.grp + K + L)
  corr.X        <- matrix(nrow = T, ncol = L * (L - 1) / 2)

  ########################################################################
  # MCMC loop
  for (t in 1:T) {

    ############################
    ## [Change 3] beta_X update — hoist invariants, use pre-computed cov
    A_tem   <- solve(t(eta_theta) %*% eta_theta) %*% t(eta_theta)
    A_tem_t <- t(A_tem)
    H_inv   <- solve(H)

    for (j in 1:M) {
      S_j     <- S_mat[j, ]
      gamma_j <- c(gamma_mat[j, ], 0)

      Cj          <- A_tem %*% betax_covmat_list[[j]] %*% A_tem_t
      Cj_inv      <- solve(Cj)
      sigma_betax <- solve(Cj_inv + H_inv)
      mu_betax    <- sigma_betax %*% (Cj_inv %*% A_tem %*% (S_j - gamma_j))
      beta_X[j, ] <- mvrnorm(n = 1, mu = mu_betax, Sigma = sigma_betax)
    }

    W <- rinvwishart(nu = M + L + 1, S = t(beta_X) %*% beta_X + diag(tau_X2))
    R <- solve(diag(sqrt(diag(W)))) %*% W %*% solve(diag(sqrt(diag(W))))
    H <- D %*% R %*% D
    corr.X[t, ] <- R[upper.tri(R)]

    #############################
    ## [Change 4a] gamma update with pre-computed Omegak_inv + vectorized nus
    tau_inv_diag <- diag(1/tau_k2, nrow = K)
    for (j in 1:M) {
      Omegak_inv_j <- Omegak_inv_list[[j]]
      nus      <- eta_theta[1:K, , drop = FALSE] %*% beta_X[j, ]
      Sigma_gk <- solve(Omegak_inv_j + tau_inv_diag)
      mu_gk    <- Sigma_gk %*% (Omegak_inv_j %*% (S_mat[j, 1:K] - nus))
      gamma_mat[j, ] <- mvrnorm(n = 1, mu = mu_gk, Sigma = Sigma_gk)
    }

    #############################
    ## [Change 4b] tau_k2 update using gamma_mat
    for (k in 1:K) {
      alpha_posterior_k <- prior_alphak[k] + M/2
      beta_posterior_k  <- prior_betak[k]  + 0.5 * sum(gamma_mat[, k]^2)
      tau_k2[k] <- 1 / rgamma(1, shape = alpha_posterior_k, rate = beta_posterior_k)
    }

    ##############################
    ## [Change 4c] eta_theta update using pre-computed eta_sub_inv + vectorized nus
    for (l in 1:L) {
      idx          <- index_lists[[l]]
      kforl        <- length(grp[[l]])
      theta_otherl <- (1:L)[-l]
      sum_precision <- matrix(0, kforl + 1, kforl + 1)
      sum_weighted  <- numeric(kforl + 1)

      for (j in 1:M) {
        bx2        <- beta_X[j, l]^2
        scaled_inv <- bx2 * eta_sub_inv[[l]][[j]]
        sum_precision <- sum_precision + scaled_inv

        if (length(theta_otherl) > 0) {
          nus <- eta_theta[idx, theta_otherl, drop = FALSE] %*% beta_X[j, theta_otherl]
        } else {
          nus <- numeric(kforl + 1)
        }
        S_j     <- S_mat[j, idx]
        gamma_j <- c(gamma_mat[j, ], 0)[idx]
        sum_weighted <- sum_weighted +
          as.numeric(scaled_inv %*% ((S_j - nus - gamma_j) / beta_X[j, l]))
      }
      Sigma_eta <- solve(sum_precision)
      mu_eta    <- Sigma_eta %*% sum_weighted
      eta_theta[idx, l] <- mvrnorm(n = 1, mu = mu_eta, Sigma = Sigma_eta)
    }

    eta_simulated[t, ] <- c(eta_theta[K+1, ],
                            unlist(lapply(1:L, function(i) eta_theta[grp[[i]], i])),
                            tau_k2)
    if (t %% 10 == 0) print(t)
  }

  ########################################################################
  # Post-processing (unchanged)
  colnames(eta_simulated) <- c(paste0('theta', 1:L),
                               unlist(lapply(1:L, function(i) paste0('theta', grp[[i]], i))),
                               paste0('h2gamma', 1:K))
  eta_simulated <- as.data.frame(eta_simulated)

  # Compute p-values for each theta_kl per latent exposure
  res.p.list <- list()
  for (l in 1:L) {
    kl_cols <- paste0("theta", grp[[l]], l)
    p_kl <- rep(0, length(kl_cols))
    ci_kl <- t(sapply(kl_cols, function(col) {
      quantile(eta_simulated[-(1:burnin), col], c(0.025, 0.975))
    }))
    for (idx2 in seq_along(kl_cols)) {
      se <- (ci_kl[idx2, 2] - ci_kl[idx2, 1]) / (2 * qnorm(0.975))
      z  <- mean(eta_simulated[-(1:burnin), kl_cols[idx2]]) / se
      tmp.p <- 2 * exp(pnorm(-abs(z), log.p = TRUE))
      p_kl[idx2] <- tmp.p
    }
    names(p_kl) <- kl_cols
    res.p.list[[l]] <- p_kl
  }

  check_sig <- sapply(1:L, function(l) any(res.p.list[[l]] < 0.3))

  if (!any(check_sig)) {
    stop("No significant thetakl found for any latent exposure. Recommend changing biomarkers.")
  }

  for (l in 1:L) {
    if (!check_sig[l]) next
    p_kl    <- res.p.list[[l]]
    idx_tmp <- which.min(p_kl)
    k_tmp   <- grp[[l]][idx_tmp]
    col_tmp <- names(p_kl)[idx_tmp]

    if (mean(eta_simulated[-(1:burnin), col_tmp]) * sign[[l]][k_tmp] < 0) {
      cols.l <- c(paste0("theta", l), paste0("theta", grp[[l]], l))
      eta_simulated[, cols.l] <- eta_simulated[, cols.l] * (-1)
    }
  }

  ci <- t(sapply(1:ncol(eta_simulated), function(x){
    quantile(eta_simulated[-(1:burnin), x], c(0.025, 0.975))
  }))
  rownames(ci) <- colnames(eta_simulated)

  corr_sub <- corr.X[-(1:burnin), ]
  if (is.vector(corr_sub)) corr_sub <- matrix(corr_sub, ncol = 1)
  cor.X <- diag(L)
  cor.X[upper.tri(cor.X)] <- colMeans(corr_sub)
  cor.X[lower.tri(cor.X)] <- t(cor.X)[lower.tri(cor.X)]

  Res_latent <- list()
  calmr_df   <- data.frame()
  for (l in 1:L) {
    if (!check_sig[l]) {
      Res_latent[[l]] <- paste0("No significant thetakl found for latent exposure X_ ", l, ". Recommend changing biomarkers.")
      next
    }
    thetal_tmp <- paste0("theta", l)

    calmr.rej <- 1
    if (ci[thetal_tmp, 1] < 0 & ci[thetal_tmp, 2] > 0) calmr.rej <- 0

    calmr.sign <- 0
    if (mean(eta_simulated[-(1:burnin), thetal_tmp]) > 0 & calmr.rej == 1) calmr.sign <-  1
    if (mean(eta_simulated[-(1:burnin), thetal_tmp]) < 0 & calmr.rej == 1) calmr.sign <- -1

    se_l <- (ci[thetal_tmp, 2] - ci[thetal_tmp, 1]) / (2 * qnorm(0.975))
    z_l  <- mean(eta_simulated[-(1:burnin), thetal_tmp]) / se_l
    p_l  <- 2 * exp(pnorm(-abs(z_l), log.p = TRUE))

    Res_latent[[l]] <- list(calmr.p = p_l, calmr.rej = calmr.rej,
                            calmr.sign = calmr.sign, thetakl.p = res.p.list[[l]])
    calmr_df <- rbind(calmr_df,
                      data.frame(calmr.p = p_l, calmr.rej = calmr.rej, calmr.sign = calmr.sign))
  }
  rownames(calmr_df) <- NULL

  Res <- list(calmr_result      = calmr_df,
              Res_latent_detail = Res_latent,
              CI                = ci,
              cor.X             = cor.X,
              mcmc.detail       = eta_simulated)
  return(Res)
}
