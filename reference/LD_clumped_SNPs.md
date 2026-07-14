# LD clump SNPs

LD clump SNPs

## Usage

``` r
LD_clumped_SNPs(
  df,
  pval_threshold,
  r2,
  kb,
  plink_path,
  reference_dir,
  outcome_name,
  create_summary_txt_file
)
```

## Arguments

- df:

  Data frame of merged biomarker + outcome SNPs (output of
  `intersection_with_outcome`).

- pval_threshold:

  P-value threshold for clumping.

- r2:

  LD r-squared clumping threshold.

- kb:

  Clumping window in kilobases.

- plink_path:

  Path to the PLINK.

- reference_dir:

  Path to reference `.bim/.bed/.fam` files.

- outcome_name:

  Character string identifying the outcome (used to exclude outcome
  p-value columns from clumping selection).

- create_summary_txt_file:

  Logical; keep the clumping summary file?

## Value

A data frame of LD-clumped SNPs.
