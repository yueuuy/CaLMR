# Select SNPs that are associated with at least two biomarkers

Select SNPs that are associated with at least two biomarkers

## Usage

``` r
significant_biomarker_snps(dfs, pval_threshold, reference_dir)
```

## Arguments

- dfs:

  A named list of biomarker GWAS data frames, each with columns `SNP`,
  `CHR`, `PVAL`, `BETA`, `SE`, `A1`, `A2`, and `NEFF`.

- pval_threshold:

  A named numeric vector of p-value thresholds, one per biomarker. Names
  must match names of `dfs`.

- reference_dir:

  Path to directory containing PLINK `.bim` files (one per chromosome,
  named `chr1.bim`, ..., `chr22.bim`).

## Value

A data frame of merged SNP-level data for SNPs passing the filter.
