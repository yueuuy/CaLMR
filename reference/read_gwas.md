# Read a GWAS Summary Statistics File

Reads GWAS summary statistics from various file formats and performs
basic standardization.

## Usage

``` r
read_gwas(path)
```

## Arguments

- path:

  File path. Supported: `.RDS`, `.RData/.rda`, `.csv`, `.tsv`, `.txt`,
  `.gz`.

## Value

A data frame with standardized column names.
