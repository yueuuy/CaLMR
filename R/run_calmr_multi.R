#' CaLMR (Multi) Pipeline: From File Paths to Multi-Exposure Causal Inference Results
#'
#' A function that reads QC'd and harmonized GWAS summary statistics from
#' file paths, selects IVSs, performs LD clumping, and runs
#' CaLMR (Multi) to jointly test the causal effects of multiple latent exposures 
#' on an outcome.
#'
#' @param biomarker_paths A named character vector of file paths to the K
#'   biomarker GWAS summary statistics files (QC'd and harmonized). Names are
#'   used as trait labels and must cover all traits referenced in \code{grp}.
#'   Supported formats: \code{.RDS}, \code{.RData}, \code{.csv}, \code{.tsv},
#'   or tab/space-delimited \code{.txt}/\code{.gz}. Each file must contain at
#'   minimum: \code{SNP}, \code{CHR}, \code{BETA}, \code{SE}, \code{PVAL},
#'   \code{A1}, \code{A2}, and \code{NEFF}.
#' @param outcome_path A single file path to the outcome GWAS summary
#'   statistics (same format requirements as biomarker files).
#' @param outcome_name A character string used as the outcome label.
#' @param traitvec A character vector of length K giving the biomarker trait
#'   names, in the order they should appear in the model. 
#'   Must match names of \code{biomarker_paths}.
#' @param grp A list of length L, where each element is a character vector of
#'   trait names belonging to that latent factor. All names must appear in
#'   \code{traitvec}.
#' @param sign A list of length L. Each element \code{sign[[l]]} is a numeric
#'   vector of length K, where the k-th element gives the pre-known sign of
#'   theta_kl for the k-th trait in \code{traitvec} (1 = positive,
#'   -1 = negative, 0 = not associated or unknown).
#' @param Corr.mat A (K+1) x (K+1) estimated correlation matrix of the biomarker 
#'   GWAS summary statistics. Row and column names must include
#'   all elements of \code{traitvec}. The outcome row/column (zeros with
#'   diagonal 1) is appended automatically inside \code{calmr_uni}.
#'   If the matrix is not provided, ld_score_dir will need to be specified
#'   for estimating the correlation matrix in the function.
#' @param pval_threshold A named numeric vector of p-value thresholds for SNP
#'   selection, one per biomarker.
#' @param reference_dir Path to the directory containing PLINK reference files
#'   (\code{chr1.bim}, ..., \code{chr22.bim}).
#' @param plink_path Path to the PLINK executable.
#' @param r2 LD clumping r-squared threshold. Default 0.001.
#' @param kb LD clumping window in kilobases. Default 500.
#' @param pvalthr_clump P-value threshold for LD clumping. Default 5e-3.
#' @param ld_score_dir Path to the directory containing LD score files for
#'   performing between-biomarker correlation. If \code{NULL}, the \code{Corr.mat}
#'   provided is used directly and no estimation will be performed.
#' @param confounder Optional character vector of SNP IDs to exclude.
#' @param T_iter Total Gibbs sampler iterations. Default 3000
#' @param burnin Burn-in iterations to discard. Default 1500.
#' @param min_snps Minimum SNPs required per latent exposure after preprocessing.
#'   Default 15. Recommend at least 15*min(K_l)
#' @param save_clump_file Logical; keep the PLINK clumping summary file?
#'   Default \code{FALSE}.
#'
#' @return A list with components (see CaLMR (Multi) function for details):
#' \describe{
#'   \item{\code{calmr_result}}{A data frame with one row per latent exposure
#'     (with at least one significant theta_kl): \code{calmr.p},
#'     \code{calmr.rej}, \code{calmr.sign}.}
#'   \item{\code{Res_latent_detail}}{A list of length L with detailed results
#'     per latent exposure.}
#'   \item{\code{CI}}{Matrix of 95\% credible intervals for all parameters.}
#'   \item{\code{cor.X}}{Estimated L x L correlation matrix among latent
#'     exposures.}
#'   \item{\code{mcmc.detail}}{Full MCMC chain (data frame).}
#'   \item{\code{calmr_info}}{Metadata: \code{n_ivs}, \code{traitvec},
#'     \code{outcome_name}, \code{grp}, \code{Corr.mat}.}
#' }
#' @export
run_calmr_multi <- function(biomarker_paths, outcome_path, outcome_name, traitvec,
                            grp, sign, Corr.mat, pval_threshold, reference_dir,
                            plink_path, r2 = 0.001, kb = 500, pvalthr_clump = 5e-3,
                            ld_score_dir = NULL, confounder = NULL,
                            T_iter = 3000, burnin = 1500, min_snps = 15,
                            save_clump_file = FALSE) {

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
      df <- data.table::fread(path, header = TRUE, showProgress = FALSE)
    }
    df <- as.data.frame(df)
    if (!"SNP" %in% colnames(df) && colnames(df)[1] != "SNP") {
      colnames(df)[1] <- "SNP"
    }
    if ("SE" %in% colnames(df)) {
      df <- df[df$SE != 0, ]
    }
    return(df)
  }

  # input checks
  L <- length(grp)
  K <- length(traitvec)

  if (L < 2) {
    stop("CaLMR(Multi) requires at least 2 factors. Use run_calmr_uni() for a single latent exposure.",
         call. = FALSE)
  }
  if (length(biomarker_paths) != K) {
    stop("Length of 'biomarker_paths' (", length(biomarker_paths), ") must equal ",
         "length of 'traitvec' (", K, ").", call. = FALSE)
  }
  if (is.null(names(biomarker_paths)) || !all(traitvec %in% names(biomarker_paths))) {
    stop("'biomarker_paths' must be a named vector with names matching 'traitvec'.",
         call. = FALSE)
  }
  # Check all grp traits exist in traitvec
  all_grp_traits <- unique(unlist(grp))
  if (!all(all_grp_traits %in% traitvec)) {
    missing <- setdiff(all_grp_traits, traitvec)
    stop("The following traits in 'grp' are not in 'traitvec': ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(sign) != L) {
    stop("'sign' must be a list of length L = ", L, ".", call. = FALSE)
  }
  for (l in seq_along(sign)) {
    if (length(sign[[l]]) != K) {
      stop("sign[[", l, "]] must have length K = ", K,
           " (one entry per trait in traitvec).", call. = FALSE)
    }
  }
  if (!all(traitvec %in% names(pval_threshold))) {
    missing <- setdiff(traitvec, names(pval_threshold))
    stop("'pval_threshold' is missing entries for: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  message("=== CaLMR (Multi) Running ===")
  message("  L = ", L, " latent exposures, K = ", K, " unique traits")
  message("  Traits: ", paste(traitvec, collapse = ", "))
  for (l in seq_along(grp)) {
    message("  Factor ", l, ": ", paste(grp[[l]], collapse = ", "))
  }

  # ---- 1. Load all GWAS data ----
  message("=== Step 1: Loading GWAS data ===")
  all_biomarkers <- lapply(biomarker_paths[traitvec], read_gwas)
  names(all_biomarkers) <- traitvec

  outcome_df <- read_gwas(outcome_path)
  message("  Loaded ", K, " biomarker GWAS and 1 outcome GWAS.")

  # Build per-factor lists of data frames from the all biomarker list
  factors <- lapply(grp, function(trait_names) all_biomarkers[trait_names])

  # ---- 2. Estimate sample overlap when needed----
  if (!is.null(ld_score_dir)) {
    message("=== Estimating correlation matrix ===")
    Corr.mat <- corr_cal( 
      gwas_list = all_biomarkers,
      traitvec = traitvec,
      outcome = outcome_name,
      ld_score_dir = ld_score_dir
    )
    message("  Finished corr.mat.")
  }

  # Validate Corr.mat
  if (!all(traitvec %in% rownames(Corr.mat))) {
    stop("Row/column names of 'Corr.mat' must include all elements of 'traitvec'.",
         call. = FALSE)
  }

  # ---- 3. IV selection per latent exposure ----
  message("=== Step 2: IV Selection per latent exposure ===")
  sample_sumtable_dfs <- list()
  
  for (f_idx in seq_along(factors)) {
    factor <- factors[[f_idx]]
    
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
  
  # ---- 4. Combine IVs ----
  message("=== Step 3: Combining and re-clumping ===")
  common_columns <- Reduce(intersect, lapply(sample_sumtable_dfs, names))
  df_list_subset <- lapply(sample_sumtable_dfs, function(df) df[, common_columns, drop = FALSE])
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
  
  # ---- 5. Run CaLMR (Multi) ----
  message("=== Step 4: Running CaLMR (Multi) ===")
  result <- calmr_multi( sumtable = combined_df, Corr.mat = Corr.mat,
                         grp = grp, L = L, K = K, sign = sign,
                         T = T_iter, burnin = burnin,
                         traitvec = traitvec, outcome = outcome_name )
  # metedata for reference
  result$calmr_info <- list( n_ivs = nrow(combined_df), 
                             traitvec = traitvec,
                             outcome_name = outcome_name,
                             grp = grp, Corr.mat = Corr.mat )
  
  message("=== Done ===")
  return(result)
}