# CaLMR (Uni) Pipeline

Performs IV selection, LD clumping, and runs CaLMR(Uni). Called by
[`run_calmr`](https://yueuuy.github.io/CaLMR/reference/run_calmr.md).
Can also be applied directly with pre-loaded data.

## Usage

``` r
run_calmr_uni(
  all_biomarkers,
  outcome_df,
  outcome_name,
  traitvec,
  sign,
  Corr.mat,
  pval_threshold,
  reference_dir,
  plink_path,
  r2,
  kb,
  pvalthr_clump,
  confounder,
  T_iter,
  burnin,
  min_snps,
  save_clump_file
)
```

## Arguments

- all_biomarkers:

  A named list of K data frames, one per biomarker.

- outcome_df:

  Data frame of outcome GWAS summary statistics.

- outcome_name:

  Character string for the outcome label.

- traitvec:

  Character vector of length K.

- sign:

  Numeric vector of length K (1, -1, or 0).

- Corr.mat:

  Correlation matrix (K x K or (K+1) x (K+1)).

- pval_threshold:

  Named numeric vector of p-value thresholds.

- reference_dir:

  Path to PLINK reference files.

- plink_path:

  Path to PLINK.

- r2:

  LD clumping r-squared threshold.

- kb:

  LD clumping window in kilobases.

- pvalthr_clump:

  P-value threshold for LD clumping.

- confounder:

  Optional character vector of SNP IDs to exclude.

- T_iter:

  Total Gibbs sampler iterations.

- burnin:

  Burn-in iterations.

- min_snps:

  Minimum total SNPs after preprocessing.

- save_clump_file:

  Logical; keep clumping summary file

## Value

A list with `calmr_result`, `CI`, `mcmc.detail`, and `calmr_info`.
Returns `NULL` if too few SNPs.
