# CaLMR

CaLMR (**C**ausal **a**nalysis of **L**atent exposures using **M**endelian **R**andomization)
is a Bayesian MR method to test the causal relationships between latent exposures and an outcome 
using GWAS summary-level association statistics of multiple traits co-regulated by the exposures. 
Both univariable (**CaLMR (Uni)**) and multivariable (**CaLMR (Multi)**) versions of CaLMR are available.

## Installation
``` R
# if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("yueuuy/CaLMR")
library(MASS)
library(LaplacesDemon)
```

## Quick Start
```r
library(CaLMR)

# ---- 1. Estimate between-biomarker correlation ----
biomarker_paths <- c(crp = "data/crp.RDS", 
                     il6 = "data/il6.RDS",
                     il8 = "data/il8.RDS")
corr_KK <- corr_cal(biomarker_paths, traitvec = c("crp", "il6", "il8"),
                    ld_score_dir = "path/to/ld_scores")

# ---- 2. Append outcome row/column ----
K <- nrow(corr_KK)
Corr.mat <- matrix(0, K + 1, K + 1) # assume a two-sample GWAS settins
Corr.mat[1:K, 1:K] <- corr_KK
Corr.mat[K + 1, K + 1] <- 1
colnames(Corr.mat) <- rownames(Corr.mat) <- c("crp", "il6", "il8", "ra")

# ---- 3. Run CaLMR ----
result <- run_calmr(
  method = "uni",
  biomarker_paths = biomarker_paths,
  outcome_path = "data/ra.RDS",
  outcome_name = "ra",
  traitvec = c("crp", "il6", "il8"),
  sign = c(1, 1, 1),
  Corr.mat = Corr.mat,
  pval_threshold = c(crp = 5e-5, il6 = 5e-5, il8 = 5e-5),
  reference_dir = "path/to/1KG/GRCh37",
  plink_path = "path/to/plink"
)

result$calmr_result
```

## Tutorial

See the full [tutorial vignette](vignettes/CaLMR_tutorial.Rmd) for a complete
walkthrough, including:

- GWAS summary statistics QC and harmonization
- Estimating the sample-overlap correlation matrix
- Running CaLMR

The tutorial demonstrates the CaLMR workflow using the Chronic Inflammation - Rheumatoid Arthritis 
analysis from the paper as an example, with CRP, IL-6, IL-8, CCL2, and TNF-α 
as the biomarkers and RA as the outcome. The package includes sample data (chromosomes 1-2) for format reference.


## Input Data Format

Each GWAS summary statistics file must contain these columns:

| Column | Description |
|--------|-------------|
| `SNP`  | SNP identifier (rsID) |
| `CHR`  | Chromosome |
| `BETA` | Effect size estimate |
| `SE`   | Standard error |
| `PVAL` | P-value |
| `A1`   | Effect allele |
| `A2`   | Other allele |
| `NEFF` | Effective sample size |

Supported file formats: `.RDS`, `.RData`, `.csv`, `.tsv`, `.txt`, `.gz`.
