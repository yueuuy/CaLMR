# Intersect Biomarker SNPs with Outcome GWAS

Intersect Biomarker SNPs with Outcome GWAS

## Usage

``` r
intersection_with_outcome(biomarker_df, outcome_df, outcome_name)
```

## Arguments

- biomarker_df:

  Data frame of biomarker SNPs (output of `significant_biomarker_snps`).

- outcome_df:

  Data frame of outcome GWAS summary statistics.

- outcome_name:

  Character string used to prefix outcome columns.

## Value

A merged data frame with SNP-level data from both sources.
