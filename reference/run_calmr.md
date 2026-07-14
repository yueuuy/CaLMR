# Run CaLMR Analysis

Complete CaLMR pipeline. Reads QC'd and harmonized GWAS summary
statistics from file paths, then directs to CaLMR (Uni) or CaLMR (Multi)
pipelines depending on the `method` argument (must specify either uni or
multi).

## Usage

``` r
run_calmr(
  method = c("uni", "multi"),
  biomarker_paths,
  outcome_path,
  outcome_name,
  traitvec,
  sign,
  Corr.mat,
  pval_threshold,
  reference_dir,
  plink_path,
  grp = NULL,
  r2 = 0.001,
  kb = 500,
  pvalthr_clump = 0.005,
  confounder = NULL,
  T_iter = 3000,
  burnin = 1500,
  min_snps = NULL,
  save_clump_file = FALSE
)
```

## Arguments

- method:

  Either `"uni"` or `"multi"`.

  `"uni"`

  :   Tests the causal effect of a single latent exposure (regulating K
      biomarkers) on the outcome.

  `"multi"`

  :   Jointly tests the causal effects of L \>= 2 latent exposures on
      the outcome. Requires `grp`.

- biomarker_paths:

  A named character vector of file paths to the K biomarker GWAS summary
  statistics files (QC'd and harmonized). Names are used as trait labels
  and must match `traitvec`. Supported formats: `.RDS`, `.RData`,
  `.csv`, `.tsv`, or tab/space- delimited `.txt`/`.gz`. Each file must
  contain at minimum: `SNP`, `CHR`, `BETA`, `SE`, `PVAL`, `A1`, `A2`,
  and `NEFF`.

- outcome_path:

  A single file path to the outcome GWAS summary statistics (same format
  requirements as biomarker files).

- outcome_name:

  A character string used as the outcome label.

- traitvec:

  A character vector of length K giving the biomarker trait names. Must
  match names of `biomarker_paths`.

- sign:

  The pre-known signs of biomarker-exposure parameters. Structure
  depends on `method` argument:

  "uni"

  :   A numeric vector of length K. Each element is 1 (positive), -1
      (negative), or 0 (unknown).

  "multi

  :   A list of length L. Each element is a numeric vector of length K,
      where the k-th value gives the sign of theta_kl for the k-th trait
      in `traitvec`.

- Corr.mat:

  A (K+1) x (K+1) estimated correlation matrix of the biomarker GWAS
  summary statistics. Row and column names must include all elements of
  `traitvec`. The outcome row/column (zeros with diagonal 1) is appended
  automatically inside `calmr_uni`.

- pval_threshold:

  A named numeric vector of p-value thresholds for IV selection, one per
  biomarker trait.

- reference_dir:

  Path to the directory containing 1000 Genomes (or equivalent)
  reference files in PLINK `.bim/.bed/.fam` format, one per chromosome
  (e.g., `chr1.bim`, ..., `chr22.bim`).

- plink_path:

  Path to the PLINK.

- grp:

  A list of length L, where each element is a character vector of trait
  names belonging to that latent factor. All names must appear in
  `traitvec`. **Required when `method = "multi"`**; ignored when
  `method = "uni"`.

- r2:

  LD clumping r-squared threshold. Default 0.001.

- kb:

  LD clumping window in kilobases. Default 500.

- pvalthr_clump:

  P-value threshold for LD clumping. Default 5e-3.

- confounder:

  Optional character vector of SNP IDs to exclude.

- T_iter:

  Total Gibbs sampler iterations. Default 3000.

- burnin:

  Burn-in iterations to discard. Default 1500.

- min_snps:

  Minimum number of SNPs required after preprocessing. Default 30 for
  uni, 15 for multi. For `"uni"`, this is the total number of IVs
  selected across all K biomarkers. For `"multi"`, this is the minimum
  required *per latent exposure* (each latent exposure must have at
  least this many SNPs after its own selection and clumping step).

- save_clump_file:

  Logical; keep the PLINK clumping summary file Default `FALSE`.

## Value

The return value of the corresponding pipeline function. See
[`run_calmr_uni`](https://yueuuy.github.io/CaLMR/reference/run_calmr_uni.md)
or
[`run_calmr_multi`](https://yueuuy.github.io/CaLMR/reference/run_calmr_multi.md)
for details.
