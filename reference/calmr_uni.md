# Univariable (single-exposure) Causal analysis of Latent exposures using Mendelian Randomization CaLMR (Uni)

Causal analysis of Latent exposures using Mendelian Randomization
(CaLMR) is an MR method that tests the causal relationships between the
outcome and the latent exposure using GWAS summary-level association
statistics. This function conducts CaLMR(Uni) test, assuming there is
one latent exposure. It is built under a two-sample MR framework and
conducts Bayesian modeling using conjugate priors and Regression with
Summary Statistics (RSS) Likelihood.

## Usage

``` r
calmr_uni(
  sumtable,
  Corr.mat,
  T = 3000,
  burnin = 1500,
  K = K,
  traitvec,
  outcome,
  sign
)
```

## Arguments

- sumtable:

  a M\*(K+1) data frame containing the GWAS summary data for K
  observable traits and the outcome.

- Corr.mat:

  a (K+1)\*(K+1) estimated correlation matrix of the GWAS summary
  statistics. The 1st-Kth variables are related to the K observable
  traits, and the last variable corresponds to the outcome. \*The order
  of the observable traits should match with the order in the 'traitvec'
  vector.

- T:

  total number of iterations for the Gibbs sampler, with default T=3000

- burnin:

  length of burn-in period, with default burnin=1500.

- K:

  number of observable traits

- traitvec:

  a vector containing the names of the observable traits. This should
  match with the column names in sumtable.

- outcome:

  the name of the outcome Y, and this should match with the column name
  in sumtable.

- sign:

  a vector of length K, where the k-th element gives the pre-known sign
  of theta_k for the k-th trait in `traitvec`. Use 1 for positive, -1
  for negative, and 0 if the sign is unknown.

## Value

- `calmr_result`: a vector contains: `calmr.p` (p-value for the causal
  effect theta), `calmr.rej` (testing for theta), and `calmr.sign` (sign
  of theta: 1, -1, or 0).

- `CI`: a matrix of posterior 95

- `mcmc.detail`: full MCMC chain for all parameters, can be used to
  check for convergence.
