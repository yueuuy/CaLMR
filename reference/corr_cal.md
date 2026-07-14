# Estimate Sample Overlap Between GWAS Biomarker Datasets Using LD Score Regression

Estimates pairwise sample overlap between K GWAS biomarker summary
statistics datasets using LD score regression intercepts. WARNING: To
use with CaLMR functions, the user should append an outcome row/column
(zeros with diagonal 1) to create a (K+1) x (K+1) matrix.

## Usage

``` r
corr_cal(biomarker_paths, traitvec, ld_score_dir, max_chi2 = 80)
```

## Arguments

- biomarker_paths:

  A named character vector of file paths to the K biomarker GWAS summary
  statistics files (QC'd and harmonized). Names are used as trait labels
  and must match `traitvec`. Supported formats: `.RDS`, `.RData`,
  `.csv`, `.tsv`, or tab/space- delimited `.txt`/`.gz`. Each file must
  contain at minimum: `SNP`, `CHR`, `BETA`, `SE`, `PVAL`, `A1`, `A2`,
  and `NEFF`.

- traitvec:

  A character vector of length K giving the names of the observable
  biomarker traits. These must match names in `gwas_list` and determine
  the ordering of rows/columns 1 to K. The order must match the order
  expected by CaLMR functions.

- ld_score_dir:

  Path to a directory containing pre-computed LD score files, one per
  chromosome. Files should be named either `<chr>.l2.ldscore.gz` or
  `chr<chr>.l2.ldscore.gz` and must contain at least the columns `SNP`
  and `L2` (or `LDSCORE`).

- max_chi2:

  Maximum allowed chi-squared statistic (\\Z^2\\) for a SNP to be
  included in the regression. SNPs exceeding this threshold are set to
  `NA`. Default is 80.

## Value

A numeric matrix of dimension \\K \times K\\.

Row and column names are set to `traitvec`.
