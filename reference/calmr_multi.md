# Multivariable (multiple exposure) Causal analysis of Latent exposures using Mendelian Randomization CaLMR (Multi)

Causal analysis of Latent exposures using Mendelian Randomization
(CaLMR) is an MR method that tests the causal relationships between the
outcome and multiple latent exposures using GWAS summary-level
association statistics. This function conducts CaLMR(Multi) test,
assuming there are multiple latent exposures. It is built under a
two-sample MR framework and conducts Bayesian modeling using conjugate
priors and Regression with Summary Statistics (RSS) Likelihood.

## Usage

``` r
calmr_multi(
  sumtable,
  Corr.mat,
  grp,
  L,
  K,
  T = 3000,
  burnin = 1500,
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

- grp:

  a list of length L with the sub-list grp\[\[l\]\] containing the names
  of the observable traits associated with l-th latent exposure.

- L:

  number of latent exposures

- K:

  number of observable traits

- T:

  total number of iterations for the Gibbs sampler, with default T=3000

- burnin:

  length of burn-in period, with default burnin=1500.

- traitvec:

  a vector containing the names of the observable traits. This should
  match with the column names in sumtable.

- outcome:

  the name of the outcome Y, and this should match with the column name
  in sumtable.

- sign:

  a list of length L. Each sub-list sign\[\[l\]\] should be a named or
  positional vector of length K, where the k-th element gives the
  pre-known sign of theta_kl for the k-th trait in `traitvec`. Use 1 for
  positive, -1 for negative, and 0 if the trait is not associated with
  the l-th latent exposure or if the sign is unknown. IMPORTANTLY,
  indexing follows the global position in `traitvec` (1 to K), not the
  local position within `grp[[l]]`.

## Value

- `calmr_result`: a data frame with one row per latent exposure with at
  least one significant corresponding theta_kl,containing three columns:
  `calmr.p` (p-value for the causal effect theta_l), `calmr.rej`
  (testing for this theta_l), and `calmr.sign` (sign of theta_l: 1, -1,
  or 0). Latent exposures with no significant theta_kl are excluded from
  this data frame.

- `Res_latent_detail`: a list of length L. For latent exposure with at
  least one significant corresponding theta_kl, each entry is a list
  containing `calmr.p`, `calmr.rej`, `calmr.sign`, and `thetakl.p` (a
  vector of p-values for each theta_kl). For latent exposures with no
  significant theta_kl, the entry is a character string recommending
  changing biomarkers.

- `CI`: a matrix of posterior 95

- `cor.X`: an L\*L estimated correlation matrix among the L latent
  exposures.

- `mcmc.detail`: full MCMC chain for all parameters, can be used to
  check for convergence.
