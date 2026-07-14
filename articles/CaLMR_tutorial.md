# CaLMR Tutorial: Test the Effect of Chronic Inflammation on Rheumatoid Arthritis

The tutorial demonstrates the CaLMR workflow using the Chronic
Inflammation - Rheumatoid Arthritis analysis from the paper as an
example, with CRP, IL-6, IL-8, CCL2, and TNF-α as the biomarkers and RA
as the outcome. The package includes sample data (chromosomes 1-2) for
format reference.

We cover four stages:

1.  Required input data format
2.  GWAS QC and harmonization (preparing your own data)
3.  Estimating the between-biomarker correlation matrix
4.  Running CaLMR

## Prerequisites

``` r
# install.packages("devtools")
devtools::install_github("your-username/CaLMR")
library(CaLMR)
```

The full pipeline additionally requires:

- [PLINK](https://www.cog-genomics.org/plink/) for LD clumping
- A 1000 Genomes reference panel in PLINK binary format
  (`chr1.bim/.bed/.fam`, …, `chr22.bim/.bed/.fam`)
- LD score files for sample-overlap estimation (e.g., from the [LDSC
  repository](https://github.com/bulik/ldsc), 1000 Genomes Phase 3
  European panel)

## Part 1: Input Data Format

CaLMR requires each GWAS summary statistics file to contain the
following columns:

| Column | Description           |
|--------|-----------------------|
| `SNP`  | SNP identifier (rsID) |
| `CHR`  | Chromosome            |
| `BETA` | Effect size estimate  |
| `SE`   | Standard error        |
| `PVAL` | P-value               |
| `A1`   | Effect allele         |
| `A2`   | Other allele          |
| `NEFF` | Effective sample size |

The package includes sample data (chromosomes 1-2 only) to illustrate
the expected format:

``` r
library(CaLMR)
data_dir <- system.file("extdata", package = "CaLMR")

# Load one biomarker as an example
crp <- readRDS(file.path(data_dir, "crp.RDS"))
head(crp)
```

    ##          SNP CHR        BP A1 A2      BETA       SE        PVAL   NEFF
    ## 10 rs1000007   2 237752054  T  C  0.003273 0.004069 0.421218076 204402
    ## 19 rs1000016   2 235690982  A  G -0.002853 0.006929 0.680601654 204402
    ## 21 rs1000017   2 235691089  C  A -0.001556 0.003659 0.670730414 204402
    ## 61 rs1000050   1 162736463  C  T  0.014668 0.005135 0.004295958 204402
    ## 66 rs1000053   2  12790328  C  T  0.003064 0.006947 0.659207288 204402
    ## 84 rs1000073   1 157255396  A  G  0.000153 0.003688 0.966896972 204402

``` r
str(crp)
```

    ## 'data.frame':    182131 obs. of  9 variables:
    ##  $ SNP : chr  "rs1000007" "rs1000016" "rs1000017" "rs1000050" ...
    ##  $ CHR : int  2 2 2 1 2 1 2 2 2 1 ...
    ##  $ BP  : int  237752054 235690982 235691089 162736463 12790328 157255396 66350062 234242347 108762878 209675578 ...
    ##  $ A1  : chr  "T" "A" "C" "C" ...
    ##  $ A2  : chr  "C" "G" "A" "T" ...
    ##  $ BETA: num  0.00327 -0.00285 -0.00156 0.01467 0.00306 ...
    ##  $ SE  : num  0.00407 0.00693 0.00366 0.00513 0.00695 ...
    ##  $ PVAL: num  0.4212 0.6806 0.6707 0.0043 0.6592 ...
    ##  $ NEFF: num  204402 204402 204402 204402 204402 ...

The sample data files included are:

- `crp.RDS`, `il6.RDS`, `il8.RDS`, `ccl2.RDS`, `tnfa.RDS` - biomarker
  GWAS
- `ra.RDS` - outcome GWAS
- `corr_biomarker.RDS` — pre-computed 5×5 between-biomarker correlation
  matrix

These are subsets (chr 1-2) for format reference only. A full analysis
requires genome-wide data.

## Part 2: GWAS QC and Harmonization

Raw GWAS files often use different column names and allele codings.
Before running CaLMR, each file must be harmonized to the format above
and aligned to a common reference genome. Below is a reference workflow
using the RA GWAS (Okada et al.) as an example.

### Loading a reference genome

``` r
library(data.table)

# Load 1000 Genomes Build 37 .bim files as the allele reference
reference_dir <- "/path/to/1KG/GRCh37"
bim_files <- file.path(reference_dir, paste0("chr", 1:22, ".bim"))
genome37 <- rbindlist(lapply(bim_files, function(f) {
  fread(f, header = FALSE,
        col.names = c("CHR", "rsid", "CM", "BP", "A1", "A2"))
}))
```

### Harmonizing a single GWAS file

Take the RA GWAS as an exmple:

``` r
# Read the raw GWAS
ra_raw <- fread("RA.GWAS.Okada.txt")

# Rename allele columns
setnames(ra_raw, c("A1", "A2"), c("effect_allele", "other_allele"))

# Merge with reference 
merged <- merge(ra_raw, genome37,
                by.x = "SNP", by.y = "rsid",
                suffixes = c(".gwas", ".ref"))

# Keep only SNPs whose alleles match the reference (either orientation)
matched <- merged[
  (effect_allele == A1 & other_allele == A2) |
  (effect_allele == A2 & other_allele == A1)
]

# Flip effect direction when alleles are swapped relative to the reference
flip_idx <- matched[effect_allele == A1 & other_allele == A2, which = TRUE]
matched[flip_idx, beta := -beta]
matched[flip_idx, c("effect_allele", "other_allele") := .(other_allele, effect_allele)]

# Select and rename to CaLMR standard format
ra_qc <- matched[, .(SNP, CHR = CHR.gwas, BP = BP,
                      A1 = effect_allele, A2 = other_allele,
                      BETA = beta, SE = se, PVAL = p, NEFF = Neff)]

saveRDS(as.data.frame(ra_qc), file = "data/ra.RDS")
```

In this analysis, the five inflammatory GWAS are processed following a
similar procedure. After this step, you should have one `.RDS` file per
trait:

    data/
      crp.RDS, il6.RDS, il8.RDS, ccl2.RDS, tnfa.RDS, ra.RDS

## Part 3: Estimating the Between-Biomarker Correlation Matrix

GWAS studies that share participants produce correlated summary
statistics even when the traits are genetically independent. CaLMR
corrects for this using a correlation matrix estimated via LD score
regression.

### Estimating with `corr_cal()`

[`corr_cal()`](https://yueuuy.github.io/CaLMR/reference/corr_cal.md)
takes file paths to the biomarker GWAS and a directory of LD score
files, and returns a K × K matrix:

``` r
traitvec <- c("crp", "il6", "il8", "ccl2", "tnfa")

biomarker_paths <- c(
  crp = "data/crp.RDS",
  il6 = "data/il6.RDS",
  il8 = "data/il8.RDS",
  ccl2 = "data/ccl2.RDS",
  tnfa = "data/tnfa.RDS"
)

corr_KK <- corr_cal(
  biomarker_paths = biomarker_paths,
  traitvec  = traitvec,
  ld_score_dir = "/path/to/ld_scores"
)
```

For this tutorial, we provide the pre-computed correlation matrix from
the paper for your reference:

``` r
data_dir <- system.file("extdata", package = "CaLMR")
traitvec <- c("crp", "il6", "il8", "ccl2", "tnfa")

corr_KK <- readRDS(file.path(data_dir, "corr_biomarker.RDS"))
print(round(corr_KK, 3))
```

    ##        crp   il6   il8  ccl2  tnfa
    ## crp  1.000 0.003 0.003 0.003 0.001
    ## il6  0.003 1.000 0.442 0.412 0.525
    ## il8  0.003 0.442 1.000 0.454 0.561
    ## ccl2 0.003 0.412 0.454 1.000 0.476
    ## tnfa 0.001 0.525 0.561 0.476 1.000

### Building the (K+1) × (K+1) matrix for CaLMR

[`run_calmr()`](https://yueuuy.github.io/CaLMR/reference/run_calmr.md)
requires a (K+1) × (K+1) matrix where the last row/column corresponds to
the outcome. Under the two-sample MR assumption (no sample overlap
between biomarker and outcome GWAS), the outcome row/column is filled
with zeros for off-diagonal entries:

``` r
K <- nrow(corr_KK)
outcome_name <- "ra"

Corr.mat <- matrix(0, nrow = K + 1, ncol = K + 1)
Corr.mat[1:K, 1:K] <- corr_KK
Corr.mat[K + 1, K + 1] <- 1
colnames(Corr.mat) <- c(traitvec, outcome_name)
rownames(Corr.mat) <- c(traitvec, outcome_name)

print(round(Corr.mat, 3))
```

    ##        crp   il6   il8  ccl2  tnfa ra
    ## crp  1.000 0.003 0.003 0.003 0.001  0
    ## il6  0.003 1.000 0.442 0.412 0.525  0
    ## il8  0.003 0.442 1.000 0.454 0.561  0
    ## ccl2 0.003 0.412 0.454 1.000 0.476  0
    ## tnfa 0.001 0.525 0.561 0.476 1.000  0
    ## ra   0.000 0.000 0.000 0.000 0.000  1

## Part 4: Running CaLMR

### Setup

Define the file paths to QC’d genome-wide GWAS data and the analysis
parameters. Both
[`corr_cal()`](https://yueuuy.github.io/CaLMR/reference/corr_cal.md) and
[`run_calmr()`](https://yueuuy.github.io/CaLMR/reference/run_calmr.md)
take the same `biomarker_paths` format - a named character vector of
file paths.

``` r
biomarker_paths <- c(
  crp = "data/crp.RDS",
  il6 = "data/il6.RDS",
  il8 = "data/il8.RDS",
  ccl2 = "data/ccl2.RDS",
  tnfa = "data/tnfa.RDS"
)
outcome_path <- "data/ra.RDS"

# P-value thresholds for IV selection (per biomarker)
pval_threshold <- c(
  crp = 5e-5,
  il6 = 1e-4,
  il8 = 1e-4,
  ccl2 = 1e-4,
  tnfa = 1e-4
)

# specify the pathways
plink_path    <- "/path/to/plink"
reference_dir <- "/path/to/1KG/GRCh37"
```

### Run the CaLMR piprline

The `sign` argument encodes prior knowledge about the direction of each
biomarker’s relationship with the latent exposure (1 = positive, -1 =
negative, 0 = unknown). Here, all five biomarkers are expected to
increase with chronic inflammation levels:

``` r
result_uni <- run_calmr(
  method = "uni",
  biomarker_paths = biomarker_paths,
  outcome_path = outcome_path,
  outcome_name = "ra",
  traitvec = traitvec,
  sign = c(1, 1, 1, 1, 1),
  Corr.mat = Corr.mat,
  pval_threshold = pval_threshold,
  reference_dir = reference_dir,
  plink_path = plink_path,
  T_iter = 10000,
  burnin = 3000
)

# Main result: p-value, rejection decision, and estimated direction
result_uni$calmr_result

# 95% credible intervals for model parameters
result_uni$CI
```

## Tips

**Choosing p-value thresholds.** The thresholds in `pval_threshold`
control how many SNPs enter the IV selection. A good starting point is
5e-6 for large GWAS.

**Minimum SNPs.** If `run_calmr` returns `NULL`, too few SNPs survived
preprocessing. Try relaxing `pval_threshold`, or reducing LD clumping
stringency (`r2`, `kb`).

**Confounder exclusion.** If you know of pleiotropic SNPs, exclude them
via the `confounder` argument.

**MCMC tuning.** The defaults (`T_iter = 3000`, `burnin = 1500`) are
reasonable for exploratory analyses. For inspection, increase to
`T_iter = 10000, burnin = 3000` and check convergence via
`result$mcmc.detail`.

## Session Info

``` r
sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] CaLMR_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.7.3          cli_3.6.6            knitr_1.51          
    ##  [4] rlang_1.3.0          xfun_0.60            otel_0.2.0          
    ##  [7] generics_0.1.4       textshaping_1.0.5    data.table_1.18.4   
    ## [10] jsonlite_2.0.0       glue_1.8.1           htmltools_0.5.9     
    ## [13] ragg_1.5.2           sass_0.4.10          rmarkdown_2.31      
    ## [16] tibble_3.3.1         evaluate_1.0.5       jquerylib_0.1.4     
    ## [19] MASS_7.3-65          fastmap_1.2.0        yaml_2.3.12         
    ## [22] lifecycle_1.0.5      compiler_4.6.1       dplyr_1.2.1         
    ## [25] fs_2.1.0             pkgconfig_2.0.3      systemfonts_1.3.2   
    ## [28] digest_0.6.39        R6_2.6.1             tidyselect_1.2.1    
    ## [31] parallel_4.6.1       pillar_1.11.1        magrittr_2.0.5      
    ## [34] bslib_0.11.0         tools_4.6.1          LaplacesDemon_16.1.8
    ## [37] pkgdown_2.2.1        cachem_1.1.0         desc_1.4.3
