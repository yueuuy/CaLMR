#' CaLMR (Multi) Pipeline
#'
#' Performs IV selection, LD clumping, and runs CaLMR(Multi).
#' Called by \code{\link{run_calmr}}. Can also be applied directly
#' with pre-loaded data.
#'
#' @param all_biomarkers A named list of K data frames, one per biomarker.
#' @param outcome_df Data frame of outcome GWAS summary statistics.
#' @param outcome_name Character string for the outcome label.
#' @param traitvec Character vector of length K.
#' @param grp List of length L; each element is a character vector of trait
#'   names for that latent factor.
#' @param sign List of length L; each element is a numeric vector of length K.
#' @param Corr.mat Correlation matrix.
#' @param pval_threshold Named numeric vector of p-value thresholds.
#' @param reference_dir Path to PLINK reference files.
#' @param plink_path Path to PLINK.
#' @param r2 LD clumping r-squared threshold.
#' @param kb LD clumping window in kilobases.
#' @param pvalthr_clump P-value threshold for LD clumping.
#' @param confounder Optional character vector of SNP IDs to exclude.
#' @param T_iter Total Gibbs sampler iterations.
#' @param burnin Burn-in iterations.
#' @param min_snps Minimum SNPs required per factor after clumping.
#' @param save_clump_file Logical; keep clumping summary file
#'
#' @return A list with \code{calmr_result}, \code{Res_latent_detail},
#'   \code{CI}, \code{cor.X}, \code{mcmc.detail}, and \code{calmr_info}.
#'   Returns \code{NULL} if too few SNPs.
#'
#' @keywords internal
run_calmr_multi <- function(all_biomarkers, outcome_df, outcome_name, traitvec,
                            grp, sign, Corr.mat, pval_threshold, reference_dir,
                            plink_path, r2, kb, pvalthr_clump, confounder,
                            T_iter, burnin, min_snps, save_clump_file) {

  K <- length(traitvec)
  L <- length(grp)

  message("=== CaLMR (Multi) ===")
  message("  L = ", L, " latent exposures, K = ", K, " traits")
  for (l in seq_along(grp)) {
    message("  Factor ", l, ": ", paste(grp[[l]], collapse = ", "))
  }

  # Build per-factor lists of data frames
  factors <- lapply(grp, function(trait_names) all_biomarkers[trait_names])

  # ---- 1. Per-factor SNP processing ----
  message("=== Step 1: IV selection per latent exposure ===")
  sample_sumtable_dfs <- list()

  for (f_idx in seq_along(factors)) {
    factor <- factors[[f_idx]]
    message("  Processing Factor ", f_idx, " (",
            paste(names(factor), collapse = ", "), ")...")

    sample_sumtable_df <- factor %>%
      significant_biomarker_snps(pval_threshold, reference_dir) %>%
      intersection_with_outcome(outcome_df, outcome_name) %>%
      LD_clumped_SNPs(pvalthr_clump, r2, kb, plink_path, reference_dir,
                      outcome_name, save_clump_file)

    if (nrow(sample_sumtable_df) < min_snps) {
      warning("Factor ", f_idx, " has only ", nrow(sample_sumtable_df),
              " SNPs after clumping (need >= ", min_snps, ").", call. = FALSE)
      return(NULL)
    }

    for (f_alt_idx in seq_along(factors)) {
      if (f_alt_idx == f_idx) next
      factoralt <- factors[[f_alt_idx]]
      named <- names(factoralt)

      for (i in seq_along(named)) {
        df_name <- named[i]
        if (any(grepl(paste0("^", df_name, "_"), colnames(sample_sumtable_df)))) next

        filtered_df <- factoralt[[i]][factoralt[[i]]$SNP %in% sample_sumtable_df$SNP, ]

        cols_to_rename <- c("BETA", "SE", "NEFF", "PVAL")
        rename_idx <- which(colnames(filtered_df) %in% cols_to_rename)
        colnames(filtered_df)[rename_idx] <-
          paste(df_name, colnames(filtered_df)[rename_idx], sep = "_")

        merge_cols <- c("SNP", paste(df_name, cols_to_rename, sep = "_"))
        merge_cols <- intersect(merge_cols, colnames(filtered_df))

        sample_sumtable_df <- merge(
          sample_sumtable_df,
          filtered_df[, merge_cols, drop = FALSE],
          by = "SNP", all.x = TRUE
        )
        sample_sumtable_df <- na.omit(sample_sumtable_df)
      }
    }

    sample_sumtable_dfs <- c(sample_sumtable_dfs, list(sample_sumtable_df))
  }

  # ---- 2. Combine factor results ----
  message("=== Step 2: Combining and re-clumping ===")
  common_columns <- Reduce(intersect, lapply(sample_sumtable_dfs, names))
  df_list_subset <- lapply(sample_sumtable_dfs,
                           function(df) df[, common_columns, drop = FALSE])
  combined_df <- do.call(rbind, df_list_subset)
  combined_df <- combined_df[!duplicated(combined_df$SNP), ]

  # confounder removal if provided
  if (!is.null(confounder)) {
    combined_df <- combined_df[!combined_df$SNP %in% confounder, ]
    message("  Removed confounder SNPs.")
  }

  # Re-clump and format
  combined_df <- combined_df %>%
    LD_clumped_SNPs(pvalthr_clump, r2, kb, plink_path, reference_dir,
                    outcome_name, save_clump_file) %>%
    convert_to_CaLMR_format(outcome_name)

  if (nrow(combined_df) < min_snps) {
    warning("Only ", nrow(combined_df), " SNPs after final processing (need >= ",
            min_snps, "). Returning NULL.", call. = FALSE)
    return(NULL)
  }
  message("  ", nrow(combined_df), " SNPs ready for CaLMR(Multi).")


  # ---- 3. Run CaLMR (Multi) ----
  message("=== Step 3: Running CaLMR (Multi) ===")
  result <- calmr_multi( sumtable = combined_df, Corr.mat = Corr.mat,
                         grp = grp, L = L, K = K, sign = sign,
                         T = T_iter, burnin = burnin,
                         traitvec = traitvec, outcome = outcome_name )

  result$calmr_info <- list( n_ivs = nrow(combined_df),
                             traitvec = traitvec,
                             outcome_name = outcome_name,
                             grp = grp, Corr.mat = Corr.mat )

  message("=== Done ===")
  return(result)
}
