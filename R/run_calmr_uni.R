#' CaLMR (Uni) Pipeline: 
#'
#' A function that reads QC'd and harmonized GWAS summary statistics from
#' file paths, performs SNP selection, LD clumping, sample-overlap correction,
#' and runs CaLMR(Uni) to test the causal effect of a latent exposure (regulate
#' K obervable biomarkers) on an outcome.
#'
#' @param biomarker_paths A named character vector of file paths to the K
#'   biomarker GWAS summary statistics files (QC'd and harmonized). Names are
#'   used as trait labels throughout the pipeline and must match \code{traitvec}.
#'   Supported formats: \code{.RDS}, \code{.RData} (first object loaded),
#'   \code{.csv}, \code{.tsv}, or tab/space-delimited \code{.txt}/\code{.gz}.
#'   Each file must contain at minimum: \code{SNP}, \code{CHR}, \code{BETA},
#'   \code{SE}, \code{PVAL}, \code{A1}, \code{A2}, and \code{NEFF}.
#' @param outcome_path A single file path to the outcome GWAS summary
#'   statistics (same format requirements as biomarker files).
#' @param outcome_name A character string used as the outcome label.
#' @param traitvec A character vector of length K giving the biomarker trait
#'   names, in the order they should appear in the correlation matrix.
#'   Must match the names of \code{biomarker_paths}.
#' @param sign A numeric vector of length K giving the pre-known signs of
#'   theta_k (1 = positive, -1 = negative, 0 = unknown).
#' @param Corr.mat A (K+1) x (K+1) estimated correlation matrix of the biomarker 
#'   GWAS summary statistics. Row and column names must include
#'   all elements of \code{traitvec}. The outcome row/column (zeros with
#'   diagonal 1) is appended automatically inside \code{calmr_uni}.
#'   If the matrix is not provided, ld_score_dir will need to be specified
#'   for estimating the correlation matrix in the function.
#' @param pval_threshold A named numeric vector of length K giving the p-value
#'   threshold for each biomarker used in IV selection. 
#'   Names must match \code{traitvec}.
#' @param reference_dir Path to the directory containing 1000 Genomes (or
#'   equivalent) reference files in PLINK \code{.bim/.bed/.fam} format, one
#'   per chromosome (e.g., \code{chr1.bim}, ..., \code{chr22.bim}).
#' @param plink_path Path to the PLINK.
#' @param r2 LD clumping r-squared threshold. Default 0.001.
#' @param kb LD clumping window in kilobases. Default 500.
#' @param pvalthr_clump P-value threshold used for LD clumping. Default 5e-3.
#' @param ld_score_dir Path to the directory containing LD score files for
#'   performing between-biomarker correlation. If \code{NULL}, the \code{Corr.mat}
#'   provided is used directly and no estimation will be performed.
#' @param confounder An optional character vector of SNP IDs to exclude
#'   (e.g., known pleiotropic SNPs). Default \code{NULL}.
#' @param T_iter Total number of Gibbs sampler iterations. Default 3000.
#' @param burnin Number of burn-in iterations to discard. Default 1500.
#' @param min_snps Minimum number of SNPs required after preprocessing to
#'   proceed with CaLMR. Default 30. Recommend at least 15*K.
#' @param save_clump_file Logical; if \code{TRUE}, keep the PLINK clumping
#'   summary file. Default \code{FALSE}.
#'
#' @return A list with components (see CaLMR (Uni) function for details):
#' \describe{
#'   \item{\code{calmr_result}}{A named vector: \code{calmr.p} (p-value),
#'     \code{calmr.rej} (reject null: 1/0), \code{calmr.sign} (direction).}
#'   \item{\code{CI}}{Matrix of 95\% credible intervals for all parameters.}
#'   \item{\code{mcmc.detail}}{Full MCMC chain (data frame).}
#'   \item{\code{calmr_info}}{A list with \code{n_ivs} (number of IVs used),
#'     \code{traitvec}, \code{outcome_name}, and \code{Corr.mat} for
#'     reproducibility.}
#' }
#' Returns \code{NULL} with a warning if fewer than \code{min_snps} SNPs
#' survive preprocessing.
#'
#' If \code{ld_score_dir} is specified, \code{corr_cal} is
#' called to estimate the K x K biomarker correlation matrix from the data,
#' overriding the supplied \code{Corr.mat}.
#'
#' @export
run_calmr_uni <- function(biomarker_paths, outcome_path, outcome_name, traitvec,
                          sign, Corr.mat, pval_threshold, reference_dir,
                          plink_path, r2 = 0.001, kb = 500,
                          pvalthr_clump = 5e-3, ld_score_dir   = NULL,
                          confounder = NULL,
                          T_iter = 3000, burnin = 1500,
                          min_snps = 30, save_clump_file = FALSE) {

  # ---- 0. Input validation ----
  K <- length(traitvec)
  if (length(biomarker_paths) != K) {
    stop("Length of 'biomarker_paths' (", length(biomarker_paths), ") must equal ",
         "length of 'traitvec' (", K, ").", call. = FALSE)
  }
  if (is.null(names(biomarker_paths)) || !all(traitvec %in% names(biomarker_paths))) {
    stop("'biomarker_paths' must be a named vector with names matching 'traitvec'.",
         call. = FALSE)
  }
  if (length(sign) != K) {
    stop("'sign' must have length K = ", K, ".", call. = FALSE)
  }
  if (!all(traitvec %in% names(pval_threshold))) {
    stop("'pval_threshold' must be a named vector with names matching 'traitvec'.",
         call. = FALSE)
  }

  # ---- Helper: read a GWAS file from path ----
  read_gwas <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "rds") {
      df <- readRDS(path)
    } else if (ext == "rdata" || ext == "rda") {
      env <- new.env(parent = emptyenv())
      load(path, envir = env)
      obj_names <- ls(env)
      if (length(obj_names) == 0) stop("No objects found in ", path, call. = FALSE)
      df <- get(obj_names[1], envir = env)
    } else if (ext == "csv") {
      df <- data.table::fread(path, sep = ",", header = TRUE, showProgress = FALSE)
    } else {
      # .txt, .tsv, .gz, or any tab/space-delimited file
      df <- data.table::fread(path, header = TRUE, showProgress = FALSE)
    }
    df <- as.data.frame(df)
    # Standardize first column to SNP if needed
    if (!"SNP" %in% colnames(df) && colnames(df)[1] != "SNP") {
      colnames(df)[1] <- "SNP"
    }
    # Basic QC: remove rows with SE == 0
    if ("SE" %in% colnames(df)) {
      df <- df[df$SE != 0, ]
    }
    return(df)
  }

  # ---- 1. Load data ----
  message("=== Loading GWAS data ===")
  all_biomarkers <- lapply(biomarker_paths[traitvec], read_gwas)
  names(all_biomarkers) <- traitvec

  outcome_df <- read_gwas(outcome_path)
  message("  Loaded ", K, " biomarker GWAS and 1 outcome GWAS.")

  # ---- 2. Estimate sample overlap ----
  if (!is.null(ld_score_dir)) {
    message("=== Estimating sample overlap from LD scores ===")
    Corr.mat <- corr_cal(
      gwas_list = all_biomarkers,
      traitvec = traitvec,
      outcome = outcome_name,
      ld_score_dir = ld_score_dir
    )
    message("  Sample overlap matrix estimated.")
  }

  # check Corr.mat dimensions
  if (!all(traitvec %in% rownames(Corr.mat))) {
    stop("Row/column names of 'Corr.mat' must include all elements of 'traitvec'.",
         call. = FALSE)
  }

  # ---- 3. SNP processing pipeline ----
  message("=== Step 1: Selecting significant biomarker SNPs ===")
  biomarker_snps <- significant_biomarker_snps(all_biomarkers, pval_threshold, reference_dir)

  if (nrow(biomarker_snps) == 0) {
    warning("No SNPs passed the biomarker p-value filter. Returning NULL.", call. = FALSE)
    return(NULL)
  }

  message("=== Step 2: Intersecting with outcome GWAS ===")
  merged_snps <- intersection_with_outcome(biomarker_snps, outcome_df, outcome_name)

  if (nrow(merged_snps) == 0) {
    warning("No SNPs in the intersection with the outcome. Returning NULL.", call. = FALSE)
    return(NULL)
  }

  message("=== Step 3: LD clumping ===")
  clumped_snps <- LD_clumped_SNPs(df = merged_snps, pval_threshold = pvalthr_clump,
                                  r2 = r2, kb = kb, plink_path = plink_path,
                                  reference_dir = reference_dir, 
                                  outcome_name = outcome_name,
                                  create_summary_txt_file = save_clump_file)

  # confounder removal if any
  if (!is.null(confounder)) {
    clumped_snps <- clumped_snps[!clumped_snps$SNP %in% confounder, ]
    message("  Removed ", length(confounder), " confounder SNPs.")
  }

  message("=== Step 4: Formatting for CaLMR ===")
  sumtable <- convert_to_CaLMR_format(clumped_snps, outcome_name)

  if (nrow(sumtable) < min_snps) {
    warning("Only ", nrow(sumtable), " SNPs after preprocessing (need >= ", min_snps,
            "). Returning NULL.", call. = FALSE)
    return(NULL)
  }
  message("  ", nrow(sumtable), " SNPs selected for running CaLMR(Uni).")

  # ---- 4. Extract the sub-correlation matrix for the traits in use ----
  #sub_Corr.mat <- extract_sub_Cov.mat(Corr.mat, traitvec)

  # ---- 4. Run CaLMR (Uni) ----
  message("=== Step 5: Running CaLMR (Uni) ===")
  result <- calmr_uni( sumtable = sumtable, Corr.mat = sub_Corr.mat,
                       T = T_iter, burnin = burnin, K = K, 
                       traitvec = traitvec, outcome = outcome_name, sign = sign )

  # metadata for reference
  result$calmr_info <- list( n_ivs = nrow(sumtable), traitvec = traitvec,
                             outcome_name = outcome_name, Corr.mat = sub_Corr.mat )

  message("=== Done ===")
  return(result)
}
