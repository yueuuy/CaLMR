#' Estimate Sample Overlap Between GWAS Biomarker Datasets Using LD Score Regression
#'
#' Estimates pairwise sample overlap between K GWAS biomarker summary statistics
#' datasets using LD score regression intercepts.
#' To use with CaLMR functions, the user should append an outcome
#' row/column (zeros with diagonal 1) to create a (K+1) x (K+1) matrix.
#'
#' @param biomarker_paths A named character vector of file paths to the K
#'   biomarker GWAS summary statistics files (QC'd and harmonized). Names are
#'   used as trait labels and must match \code{traitvec}. Supported formats:
#'   \code{.RDS}, \code{.RData}, \code{.csv}, \code{.tsv}, or tab/space-
#'   delimited \code{.txt}/\code{.gz}. Each file must contain at minimum:
#'   \code{SNP}, \code{CHR}, \code{BETA}, \code{SE}, \code{PVAL}, \code{A1},
#'   \code{A2}, and \code{NEFF}.
#' @param traitvec A character vector of length K giving the names of the
#'   observable biomarker traits. These must match names in \code{gwas_list}
#'   and determine the ordering of rows/columns 1 to K. The order must match
#'   the order expected by CaLMR functions.
#' @param ld_score_dir Path to a directory containing pre-computed LD score
#'   files, one per chromosome. Files should be named either
#'   \code{<chr>.l2.ldscore.gz} or \code{chr<chr>.l2.ldscore.gz} and must
#'   contain at least the columns \code{SNP} and \code{L2} (or \code{LDSCORE}).
#' @param max_chi2 Maximum allowed chi-squared statistic (\eqn{Z^2}) for a SNP
#'   to be included in the regression. SNPs exceeding this threshold are set to
#'   \code{NA}. Default is 80.
#'
#' @return A numeric matrix of dimension \eqn{K \times K}.
#'   \describe{
#'      Biomarker traits in the order of \code{traitvec}. Off-diagonal entries
#'      are the estimated LDSC regression intercepts representing sample overlap.
#'   }
#'   Row and column names are set to \code{traitvec}.


corr_cal <- function(biomarker_paths, traitvec, ld_score_dir, max_chi2 = 80) {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("The 'data.table' package is required. Please install it.", call. = FALSE)
    }

    K <- length(traitvec)

    # check inputs
    if (length(biomarker_paths) != K) {
     stop("Length of 'biomarker_paths' (", length(biomarker_paths), ") must equal ",
         "length of 'traitvec' (", K, ").", call. = FALSE)
    }
    if (is.null(names(biomarker_paths)) || !all(traitvec %in% names(biomarker_paths))) {
      stop("'biomarker_paths' must be a named vector with names matching 'traitvec'.",
          call. = FALSE)
    }

    # load gwas data
    message("Loading biomarker GWAS data...")
    gwas_list <- lapply(biomarker_paths[traitvec], read_gwas)
    names(gwas_list) <- traitvec


    # load LD scores
    load_ld_scores <- function(dir_path, chromosomes = 1:22) {
      message("-> Loading LD scores from: ", dir_path)
      ld_scores_list <- list()
      for (chr in chromosomes) {
        ld_file <- file.path(dir_path, paste0(chr, ".l2.ldscore.gz"))
        if (!file.exists(ld_file)) {
          ld_file <- file.path(dir_path, paste0("chr", chr, ".l2.ldscore.gz"))
          if (!file.exists(ld_file)) {
            warning("LD score file for chromosome ", chr, " not found. Skipping.",
                    call. = FALSE)
            next
          }
        }
        ld_data <- data.table::fread(ld_file, showProgress = FALSE)
        if ("L2" %in% names(ld_data)) {
          data.table::setnames(ld_data, "L2", "LDSCORE")
        }
        ld_scores_list[[chr]] <- ld_data[, c("SNP", "LDSCORE")]
      }
      if (length(ld_scores_list) == 0) {
        stop("No LD score files were successfully loaded.", call. = FALSE)
      }
      combined_ld <- data.table::rbindlist(ld_scores_list)
      message("-> Successfully loaded ", format(nrow(combined_ld), big.mark = ","),
              " SNPs from LD score files.")
      return(combined_ld)
    }

    ld_scores <- load_ld_scores(ld_score_dir)

    # ---- Step 1: Common SNPs ----
    message("Step 1: Finding common SNPs across all GWAS datasets and LD scores...")
    snp_lists <- lapply(gwas_list, `[[`, "SNP")
    common_snps <- Reduce(intersect, snp_lists)
    common_snps <- intersect(common_snps, ld_scores$SNP)

    # --- FIX 1: Filter LD scores to the common set before further QC ---
    ld_scores_subset <- ld_scores[match(common_snps, ld_scores$SNP), ]
    # --- FIX 2: Add filter for extreme LD scores, matching the original script ---
    ld_qc_pass <- ld_scores_subset$LDSCORE <
      stats::quantile(ld_scores_subset$LDSCORE, 0.99, na.rm = TRUE)
    common_snps_qc <- ld_scores_subset$SNP[ld_qc_pass]
    ld_values_qc   <- ld_scores_subset$LDSCORE[ld_qc_pass]

    message("-> Found ", format(length(common_snps_qc), big.mark = ","),
            " common SNPs after QC for analysis.")

    if (length(common_snps_qc) < 50000) {
      warning("Fewer than 50,000 common SNPs found. Results may be less reliable.",
              call. = FALSE)
    }

    # ---- Step 2: Z-score matrix (K biomarkers only) ----
    message("Step 2: Preparing Z-score matrix...")
    z_matrix <- matrix(NA, nrow = length(common_snps_qc), ncol = K)
    colnames(z_matrix) <- traitvec

    for (i in seq_along(gwas_list)) {
      gwas_data <- as.data.frame(gwas_list[[i]])
      gwas_subset <- gwas_data[match(common_snps_qc, gwas_data$SNP), ]
      z_scores <- gwas_subset$BETA / gwas_subset$SE
      z_scores[z_scores^2 > max_chi2] <- NA
      z_matrix[, i] <- z_scores
    }

    # ---- Step 3: Pairwise regressions (K x K biomarker block) ----
    message("Step 3: Running pairwise regressions to estimate sample overlap...")
    overlap_KK <- matrix(NA, nrow = K, ncol = K)

    pb <- utils::txtProgressBar(min = 0, max = K * (K + 1) / 2, style = 3)
    progress_counter <- 0

    for (i in 1:K) {
      for (j in i:K) {
        if (i == j) {
          overlap_KK[i, j] <- 1.0
        } else {
          z_product <- z_matrix[, i] * z_matrix[, j]
          valid_indices <- !is.na(z_product) & !is.na(ld_values_qc)
          # --- FIX 3: Define the weights for the regression ---
          # This is the standard weighting scheme used in LDSC.
          weights <- 1 / (2 * ld_values_qc[valid_indices] + 1)
          # --- FIX 4: Add the 'weights' argument to the lm() call ---
          lm_fit <- stats::lm(z_product[valid_indices] ~ ld_values_qc[valid_indices],
                              weights = weights)
          intercept <- stats::coef(lm_fit)[1]

          overlap_KK[i, j] <- intercept
          overlap_KK[j, i] <- intercept
        }
        progress_counter <- progress_counter + 1
        utils::setTxtProgressBar(pb, progress_counter)
      }
    }
    close(pb)

    # name the row / column names
    rownames(overlap_KK) <- traitvec
    colnames(overlap_KK) <- traitvec

    message("\nDone. Returning the between-biomarker correlation matrix.")


    return(overlap_KK)
  }
