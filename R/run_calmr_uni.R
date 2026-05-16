#' CaLMR (Uni) Pipeline
#'
#' Performs IV selection, LD clumping, and runs CaLMR(Uni). Called by
#' \code{\link{run_calmr}}. Can also be applied directly with pre-loaded data.
#'
#' @param all_biomarkers A named list of K data frames, one per biomarker.
#' @param outcome_df Data frame of outcome GWAS summary statistics.
#' @param outcome_name Character string for the outcome label.
#' @param traitvec Character vector of length K.
#' @param sign Numeric vector of length K (1, -1, or 0).
#' @param Corr.mat Correlation matrix (K x K or (K+1) x (K+1)).
#' @param pval_threshold Named numeric vector of p-value thresholds.
#' @param reference_dir Path to PLINK reference files.
#' @param plink_path Path to PLINK.
#' @param r2 LD clumping r-squared threshold.
#' @param kb LD clumping window in kilobases.
#' @param pvalthr_clump P-value threshold for LD clumping.
#' @param confounder Optional character vector of SNP IDs to exclude.
#' @param T_iter Total Gibbs sampler iterations.
#' @param burnin Burn-in iterations.
#' @param min_snps Minimum total SNPs after preprocessing.
#' @param save_clump_file Logical; keep clumping summary file
#'
#' @return A list with \code{calmr_result}, \code{CI}, \code{mcmc.detail},
#'   and \code{calmr_info}. Returns \code{NULL} if too few SNPs.
#'
#' @keywords internal
run_calmr_uni <- function(all_biomarkers, outcome_df, outcome_name, traitvec,
                          sign, Corr.mat, pval_threshold, reference_dir,
                          plink_path, r2, kb, pvalthr_clump, confounder,
                          T_iter, burnin, min_snps, save_clump_file) {

  K <- length(traitvec)

  # ---- 1. SNP processing ----
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

  # ---- 2. Run CaLMR (Uni) ----
  message("=== Step 5: Running CaLMR (Uni) ===")
  result <- calmr_uni( sumtable = sumtable, Corr.mat = sub_Corr.mat,
                       T = T_iter, burnin = burnin, K = K,
                       traitvec = traitvec, outcome = outcome_name, sign = sign )

  result$calmr_info <- list( n_ivs = nrow(sumtable), traitvec = traitvec,
                             outcome_name = outcome_name, Corr.mat = sub_Corr.mat )

  message("=== Done ===")
  return(result)
}
