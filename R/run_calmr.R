#' Run CaLMR Analysis
#'
#' Complete CaLMR pipeline.
#' Reads QC'd and harmonized GWAS summary statistics from file paths,
#' then directs to CaLMR (Uni) or CaLMR (Multi) pipelines depending on
#' the \code{method} argument (must specify either uni or multi).
#'
#' @param method Either \code{"uni"} or \code{"multi"}.
#'   \describe{
#'     \item{\code{"uni"}}{Tests the causal effect of a single latent exposure
#'       (regulating K biomarkers) on the outcome.}
#'     \item{\code{"multi"}}{Jointly tests the causal effects of L >= 2 latent
#'       exposures on the outcome. Requires \code{grp}.}
#'   }
#' @param biomarker_paths A named character vector of file paths to the K
#'   biomarker GWAS summary statistics files (QC'd and harmonized). Names are
#'   used as trait labels and must match \code{traitvec}. Supported formats:
#'   \code{.RDS}, \code{.RData}, \code{.csv}, \code{.tsv}, or tab/space-
#'   delimited \code{.txt}/\code{.gz}. Each file must contain at minimum:
#'   \code{SNP}, \code{CHR}, \code{BETA}, \code{SE}, \code{PVAL}, \code{A1},
#'   \code{A2}, and \code{NEFF}.
#' @param outcome_path A single file path to the outcome GWAS summary
#'   statistics (same format requirements as biomarker files).
#' @param outcome_name A character string used as the outcome label.
#' @param traitvec A character vector of length K giving the biomarker trait
#'   names. Must match names of \code{biomarker_paths}.
#' @param sign The pre-known signs of biomarker-exposure parameters.
#' Structure depends on \code{method} argument:
#'   \describe{
#'     \item{"uni"}{A numeric vector of length K. Each element is 1 (positive),
#'       -1 (negative), or 0 (unknown).}
#'     \item{"multi}{A list of length L. Each element is a numeric vector of
#'       length K, where the k-th value gives the sign of theta_kl for the
#'       k-th trait in \code{traitvec}.}
#'   }
#' @param Corr.mat A (K+1) x (K+1) estimated correlation matrix of the biomarker
#'   GWAS summary statistics. Row and column names must include
#'   all elements of \code{traitvec}. The outcome row/column (zeros with
#'   diagonal 1) is appended automatically inside \code{calmr_uni}.
#'   If the matrix is not provided, ld_score_dir will need to be specified
#'   for estimating the correlation matrix in the function.
#' @param pval_threshold A named numeric vector of p-value thresholds for IV
#'   selection, one per biomarker trait.
#' @param reference_dir Path to the directory containing 1000 Genomes (or
#'   equivalent) reference files in PLINK \code{.bim/.bed/.fam} format, one
#'   per chromosome (e.g., \code{chr1.bim}, ..., \code{chr22.bim}).
#' @param plink_path Path to the PLINK.
#' @param grp A list of length L, where each element is a character vector of
#'   trait names belonging to that latent factor. All names must appear in
#'   \code{traitvec}.
#'   \strong{Required when \code{method = "multi"}};
#'            ignored when \code{method = "uni"}.
#' @param r2 LD clumping r-squared threshold. Default 0.001.
#' @param kb LD clumping window in kilobases. Default 500.
#' @param pvalthr_clump P-value threshold for LD clumping. Default 5e-3.
#' @param ld_score_dir Path to the directory containing LD score files for
#'   performing between-biomarker correlation. If \code{NULL}, the \code{Corr.mat}
#'   provided is used directly and no estimation will be performed.
#' @param confounder Optional character vector of SNP IDs to exclude.
#' @param T_iter Total Gibbs sampler iterations. Default 3000.
#' @param burnin Burn-in iterations to discard. Default 1500.
#' @param min_snps Minimum number of SNPs required after preprocessing.
#'   Default 30 for uni, 15 for multi.
#'   For \code{"uni"}, this is the total number of IVs selected across all K
#'   biomarkers. For \code{"multi"}, this is the minimum required
#'   \emph{per latent exposure} (each latent exposure must have at least this many
#'   SNPs after its own selection and clumping step).
#' @param save_clump_file Logical; keep the PLINK clumping summary file
#'   Default \code{FALSE}.
#'
#' @return The return value of the corresponding pipeline function. See
#'   \code{\link{run_calmr_uni}} or \code{\link{run_calmr_multi}} for details.
#'
#' @export
run_calmr <- function(method = c("uni", "multi"),
                      biomarker_paths, outcome_path, outcome_name, traitvec,
                      sign, Corr.mat = NULL,
                      pval_threshold, reference_dir, plink_path,
                      grp = NULL, r2 = 0.001, kb = 500, pvalthr_clump = 5e-3,
                      ld_score_dir = NULL, confounder = NULL,
                      T_iter = 3000, burnin = 1500, min_snps = NULL,
                      save_clump_file = FALSE) {

  method <- match.arg(method)
  K <- length(traitvec)

  # ---- overall check ----
  if (length(biomarker_paths) != K) {
    stop("Length of 'biomarker_paths' (", length(biomarker_paths), ") must equal ",
         "length of 'traitvec' (", K, ").", call. = FALSE)
  }
  if (is.null(names(biomarker_paths)) || !all(traitvec %in% names(biomarker_paths))) {
    stop("'biomarker_paths' must be a named vector with names matching 'traitvec'.",
         call. = FALSE)
  }
  if (!all(traitvec %in% names(pval_threshold))) {
    stop("'pval_threshold' must be a named vector with names matching 'traitvec'.",
         call. = FALSE)
  }
  if (is.null(Corr.mat) && is.null(ld_score_dir)) {
    stop("Either 'Corr.mat' or 'ld_score_dir' must be provided.", call. = FALSE)
  }
  if (!is.null(Corr.mat)) {
    Corr.mat <- as.matrix(Corr.mat)
    if (!all(dim(Corr.mat) == c(K + 1, K + 1))) {
      stop( "'Corr.mat' must be a square matrix with dimensions (K + 1) x (K + 1),
            where K is the number of biomarkers. ", call. = FALSE)
    }
  }

  # ---- method-specific check ----
  if (method == "uni") {
    if (!is.numeric(sign) || !is.vector(sign) || length(sign) != K) {
      stop("For method = 'uni', 'sign' must be a numeric vector of length K = ", K, ".",
           call. = FALSE)
    }
    if (is.null(min_snps)) min_snps <- 30
  } else {
    # multi
    if (is.null(grp)) {
      stop("For method = 'multi', 'grp' must be provided.", call. = FALSE)
    }
    L <- length(grp)
    if (L < 2) {
      stop("CaLMR(Multi) requires at least 2 latent exposures. Use method = 'uni' for a single latent exposure.",
           call. = FALSE)
    }

    if (!is.list(sign) || length(sign) != L) {
      stop("For method = 'multi', 'sign' must be a list of length L = ", L, ".",
           call. = FALSE)
    }
    for (l in seq_along(sign)) {
      if (length(sign[[l]]) != K) {
        stop("sign[[", l, "]] must have length K = ", K, ".", call. = FALSE)
      }
    }
    if (is.null(min_snps)) min_snps <- 15
  }

  # ---- read gwas ----
  message("=== Loading GWAS data ===")
  all_biomarkers <- lapply(biomarker_paths[traitvec], read_gwas)
  names(all_biomarkers) <- traitvec

  outcome_df <- read_gwas(outcome_path)
  message("  Loaded GWAS.")

  # ---- estimate correlation matrix  ----
  if (!is.null(ld_score_dir)) {
    message("=== Estimating correlation matrix from LD scores ===")
    Corr.mat <- corr_cal(
      gwas_list = all_biomarkers,
      traitvec = traitvec,
      outcome = outcome_name,
      ld_score_dir = ld_score_dir
    )
    message("  Correlation matrix estimated.")
  }

  if (!all(traitvec %in% rownames(Corr.mat))) {
    stop("Names of 'Corr.mat' must include all elements of 'traitvec'.",
         call. = FALSE)
  }

  # ---- Route to method ----
  if (method == "uni") {
    result <- run_calmr_uni(
      all_biomarkers = all_biomarkers,
      outcome_df = outcome_df,
      outcome_name = outcome_name,
      traitvec = traitvec,
      sign = sign,
      Corr.mat = Corr.mat,
      pval_threshold = pval_threshold,
      reference_dir = reference_dir,
      plink_path = plink_path,
      r2 = r2, kb = kb,
      pvalthr_clump = pvalthr_clump,
      confounder = confounder,
      T_iter = T_iter, burnin = burnin,
      min_snps = min_snps,
      save_clump_file = save_clump_file
    )
  } else {
    result <- run_calmr_multi(
      all_biomarkers = all_biomarkers,
      outcome_df = outcome_df,
      outcome_name = outcome_name,
      traitvec = traitvec,
      grp = grp, sign = sign,
      Corr.mat = Corr.mat,
      pval_threshold = pval_threshold,
      reference_dir = reference_dir,
      plink_path = plink_path,
      r2 = r2, kb = kb,
      pvalthr_clump = pvalthr_clump,
      confounder = confounder,
      T_iter = T_iter, burnin = burnin,
      min_snps = min_snps,
      save_clump_file = save_clump_file
    )
  }

  return(result)
}

