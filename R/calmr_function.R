# SNP Processing Functions -------------------------------------------------
#' Package imports
#'
#' @import dplyr
#' @importFrom data.table fread
#' @importFrom tools file_ext
#' @importFrom magrittr %>%
#' @keywords internal
"_PACKAGE"
# STEP 1: Function to find the set of SNPs with p-values that fall below a threshold for at least 2 biomarkers,
# also excludes SNPs that are not present in Genome Build 37 and imputes BP
# Note: reference_dir is the directory where you store genome build 37 data

#' Select SNPs that are associated with at least two biomarkers
#'
#' @param dfs A named list of biomarker GWAS data frames, each with columns
#'   \code{SNP}, \code{CHR}, \code{PVAL}, \code{BETA}, \code{SE}, \code{A1},
#'   \code{A2}, and \code{NEFF}.
#' @param pval_threshold A named numeric vector of p-value thresholds, one per
#'   biomarker. Names must match names of \code{dfs}.
#' @param reference_dir Path to directory containing PLINK \code{.bim} files
#'   (one per chromosome, named \code{chr1.bim}, ..., \code{chr22.bim}).
#'
#' @return A data frame of merged SNP-level data for SNPs passing the filter.
#' @keywords internal
significant_biomarker_snps <- function(dfs, pval_threshold, reference_dir) {
  # Load SNPs and BPs from reference files into a hash table
  load_snps_from_bim_files <- function(reference_dir) {
    snp_hash <- new.env(hash = TRUE, parent = emptyenv())
    for (chr in 1:22) {
      bim_file <- file.path(reference_dir, paste0("chr", chr, ".bim"))
      if (file.exists(bim_file)) {
        bim_data <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
        snps <- bim_data$V2  # Assuming the SNP column is the second column in .bim files
        bps <- bim_data$V4  # Assuming the BP column is the fourth column in .bim files
        for (i in seq_along(snps)) {
          snp_hash[[snps[i]]] <- bps[i]
        }
      } else {
        warning("BIM file for chromosome ", chr, " does not exist in the reference directory.")
      }
    }
    return(snp_hash)
  }

  snp_hash <- load_snps_from_bim_files(reference_dir)

  # Prepare and filter data frames
  prepared_dfs <- lapply(names(dfs), function(name) {
    df <- dfs[[name]]
    # Standardize SNP identifiers if necessary
    if (any(grepl("^chr", df$SNP))) {
      df$SNP <- gsub("^chr[0-9XY]+_", "", df$SNP)
    }
    # Filter for SNPs starting with 'rs'
    df <- df[grep("^rs", df$SNP), ]
    # Check for necessary columns
    if (!all(c("SNP", "CHR", "PVAL") %in% colnames(df))) {
      stop("Data frame ", name, " is missing one of the required columns: SNP, CHR, PVAL.")
    }
    # Filter to include only relevant rows based on PVAL and mark the source
    df <- df[df$PVAL < pval_threshold[name], c("SNP", "CHR")]
    if (nrow(df) > 0) {
      df$source = name
    }
    df
  })

  # Combine all SNP entries from all data frames
  all_snps <- do.call(rbind, prepared_dfs)

  # Ensure there are entries to process
  if (nrow(all_snps) == 0) {
    return(data.frame())  # Return empty data frame if no SNPs qualify
  }

  # Aggregate to find SNPs present with low PVAL in at least two data frames
  snp_counts <- aggregate(source ~ SNP + CHR, data = all_snps, FUN = function(x) paste(unique(x), collapse = ", "))
  snp_counts <- snp_counts[sapply(snp_counts$source, function(x) length(unlist(strsplit(x, ",")))) >= 2, ]

  # Merge data frames including all columns, tagging columns by dataframe
  merged_df <- NULL
  for (name in names(dfs)) {
    df <- dfs[[name]]
    # Filter to include only common SNPs
    df <- df[df$SNP %in% snp_counts$SNP & df$CHR %in% snp_counts$CHR, ]
    # Rename columns to include dataframe name as prefix, skip SNP and CHR
    colnames(df)[!colnames(df) %in% c("SNP", "CHR")] <- paste(name, colnames(df)[!colnames(df) %in% c("SNP", "CHR")], sep = "_")
    # Merge or bind as appropriate
    if (is.null(merged_df)) {
      merged_df <- df
    } else {
      merged_df <- merge(merged_df, df, by = c("SNP", "CHR"), all = TRUE)
    }
  }

  # Adding a column for data frames meeting the PVAL threshold
  if (!is.null(merged_df) && nrow(merged_df) > 0) {
    merged_df$Biomarkers_Meeting_PVAL <- snp_counts$source[match(paste(merged_df$SNP, merged_df$CHR), paste(snp_counts$SNP, snp_counts$CHR))]
  }

  # Ensure at least two valid PVAL entries per row in the merged dataframe
  if (!is.null(merged_df)) {
    pval_cols <- grep("_PVAL$", colnames(merged_df), value = TRUE)
    merged_df <- merged_df[rowSums(!is.na(merged_df[, pval_cols])) >= 2, ]
  }

  # Extract columns that have SNP, CHR, BP, PVAL, BETA, SE, A1, or A2 in their names
  cols_to_keep <- grep("SNP|CHR|BP|PVAL|NEFF|BETA|SE|A1|A2", colnames(merged_df), value = TRUE)
  cols_to_remove <- grep("_BP", colnames(merged_df), value = TRUE)

  # Step 3: Subset the data frame to keep only the desired columns
  cols_to_keep <- setdiff(cols_to_keep, cols_to_remove)
  merged_df <- merged_df[, cols_to_keep, drop = FALSE]

  # Filter SNPs based on presence in the reference files and add BP information
  if (!is.null(merged_df) && nrow(merged_df) > 0) {
    merged_df <- merged_df[sapply(merged_df$SNP, function(snp) {
      exists(snp, envir = snp_hash, inherits = FALSE)
    }), ]
    # Add BP column
    merged_df$BP <- sapply(merged_df$SNP, function(snp) {
      if (exists(snp, envir = snp_hash, inherits = FALSE)) {
        return(snp_hash[[snp]])
      } else {
        return(NA)
      }
    })
  }

  if (nrow(merged_df) == 0) {
    warning('No SNPs met the specified p-value threshold for 2 or more biomarkers. Consider raising the threshold to continue with your analysis.', call. = FALSE)
  }

  merged_df <- merged_df %>%
    dplyr::select(1:2, BP, everything())


  return(merged_df)
}




# STEP 2: Function to find the intersections of SNPs between the thresholded biomarker GWAS and the outcome GWAS
# outcome_name is string indicating the name of the outcome
#' Intersect Biomarker SNPs with Outcome GWAS
#'
#' @param biomarker_df Data frame of biomarker SNPs (output of
#'   \code{significant_biomarker_snps}).
#' @param outcome_df Data frame of outcome GWAS summary statistics.
#' @param outcome_name Character string used to prefix outcome columns.
#'
#' @return A merged data frame with SNP-level data from both sources.
#' @keywords internal
intersection_with_outcome <- function(biomarker_df, outcome_df, outcome_name) {
  # Ensure necessary columns are present
  if (nrow(biomarker_df) == 0) {
    stop('Biomarker dataframe contains no SNPs.')
  }

  if (nrow(outcome_df) == 0) {
    stop('Outcome dataframe contains no SNPs.')
  }


  # Extract columns that have SNP, CHR, BP, PVAL, BETA, SE, A1, or A2 in their names
  cols_to_keep <- grep("PVAL|SNP|CHR|BETA|SE|A1|A2|NEFF", colnames(outcome_df), value = TRUE)

  outcome_df <- outcome_df[, cols_to_keep, drop = FALSE]


  if (!all(c("SNP", "CHR") %in% colnames(biomarker_df)) || !all(c("SNP", "CHR") %in% colnames(outcome_df))) {
    stop("Both data frames must contain SNP and CHR columns.")
  }

  # Merge data frames based on SNP and CHR to find intersection
  merged_df <- merge(biomarker_df, outcome_df, by = c("SNP", "CHR"), all = FALSE)

  # If you need to rename the columns from solo_df to indicate their source
  outcome_cols <- setdiff(colnames(outcome_df), c("SNP", "CHR"))
  colnames(merged_df)[colnames(merged_df) %in% outcome_cols] <- paste(outcome_name, colnames(merged_df)[colnames(merged_df) %in% outcome_cols], sep = "_")


  if (nrow(merged_df) == 0) {
    warning('No SNPs were found in the intersection with the outcome. Consider raising the p-value threshold in the previous function to continue with your analysis.', call. = F)
  }


  return(merged_df)
}




# STEP 3: Function to perform LD clumping
# create_summary_txt_file is a boolean when set to true will produce a .txt file in the current wd containing a summary of the clumped SNPs
# plink_path is directory where software for plink is stored
#' LD clump SNPs
#' @param df Data frame of merged biomarker + outcome SNPs (output of
#'   \code{intersection_with_outcome}).
#' @param pval_threshold P-value threshold for clumping.
#' @param r2 LD r-squared clumping threshold.
#' @param kb Clumping window in kilobases.
#' @param plink_path Path to the PLINK.
#' @param reference_dir Path to reference \code{.bim/.bed/.fam} files.
#' @param outcome_name Character string identifying the outcome (used to
#'   exclude outcome p-value columns from clumping selection).
#' @param create_summary_txt_file Logical; keep the clumping summary file?
#'
#' @return A data frame of LD-clumped SNPs.
#' @keywords internal
#'
LD_clumped_SNPs <- function(df, pval_threshold, r2, kb, plink_path, reference_dir, outcome_name, create_summary_txt_file) {

  if (nrow(df) == 0) {
    stop('No SNPs were found in the intersection with the outcome. Consider raising the p-value threshold in the first function to continue with your analysis.', call. = F)
  }

  # Extract necessary columns based on pattern and exclude outcome_ columns
  pval_cols <- grep(paste0("^(?!", outcome_name, "_).*_PVAL$"), colnames(df), value = TRUE, perl = TRUE)

  # Exclude the Biomarkers_Meeting_PVAL column
  pval_cols <- setdiff(pval_cols, "Biomarkers_Meeting_PVAL")

  # Ensure pval columns are numeric
  df <- df %>%
    dplyr::mutate(across(all_of(pval_cols), as.numeric))

  # Filter for p-values below the threshold
  pval_filter_expr <- paste(pval_cols, "<", pval_threshold, collapse = " | ")
  filtered_snps <- df %>%
    dplyr::filter(eval(parse(text = pval_filter_expr)))

  # Select the largest p-value below the threshold for each SNP
  filtered_snps <- filtered_snps %>%
    rowwise() %>%
    dplyr::mutate(max_pval_below_threshold = max(c_across(all_of(pval_cols))[c_across(all_of(pval_cols)) < pval_threshold], na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::filter(!is.na(max_pval_below_threshold))

  # Ensure that the correct columns are selected
  bp_col <- grep(paste0("^(?!", outcome_name, "_).*_BP$"), colnames(filtered_snps), value = TRUE, perl = TRUE)


  # Select only necessary columns for clumping
  snps_for_clumping <- filtered_snps %>%
    dplyr::select(SNP, CHR, BP, max_pval_below_threshold)

  # Rename columns to match PLINK requirements
  colnames(snps_for_clumping) <- c("SNP", "CHR", "BP", "P")

  # Save to a file for PLINK
  write.table(snps_for_clumping, "snps_for_clumping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

  # Define the PLINK clumping command for a single chromosome
  run_plink_clumping <- function(chr) {
    cmd <- sprintf("%s --bfile %s/chr%d --clump snps_for_clumping.txt --clump-p1 %s --clump-r2 %s --clump-kb %s --out clumped_snps_chr%d",
                   plink_path, reference_dir, chr, pval_threshold, r2, kb, chr)
    system(cmd)
  }

  # Loop over chromosomes and run PLINK clumping
  for (chr in 1:22) {
    run_plink_clumping(chr)
  }

  # Combine all clumped files into a single file
  clumped_files <- list.files(pattern = "clumped_snps_chr\\d+\\.clumped")
  combined_clumped <- do.call(rbind, lapply(clumped_files, read.table, header = TRUE))
  write.table(combined_clumped, "summary_SNPs_clumped.txt", sep = "\t", row.names = FALSE, quote = FALSE)

  # Delete temporary files
  temp_files <- c(list.files(pattern = "clumped_snps_chr\\d+\\.log"),
                  list.files(pattern = "clumped_snps_chr\\d+\\.nosex"),
                  clumped_files)
  file.remove(temp_files)
  file.remove('snps_for_clumping.txt')


  # Load the combined clumped results
  clumped_snps <- read.table("summary_SNPs_clumped.txt", header = TRUE, sep = "\t")

  # Extract a subset of rows from the original dataframe based on the rsids in clumped results
  subset_df <- df %>%
    dplyr::filter(SNP %in% clumped_snps$SNP)

  if (!create_summary_txt_file) {
    file.remove('summary_SNPs_clumped.txt')
  }


  # Return the subset dataframe
  return(as.data.frame(subset_df))
}




# Dataframe Formatting Function for CaLMR
#' Convert merged GWAS data to CaLMR input format
#'
#' @param df Data frame of merged SNP-level data (output of
#'   \code{LD_clumped_SNPs}).
#' @param outcome_trait Character string identifying the outcome.
#' @keywords internal
convert_to_CaLMR_format <- function(df, outcome_trait) {
  # Extract columns with *_BETA, *_SE, and *_NEFF
  beta_cols <- grep("_BETA$", names(df), value = TRUE)
  se_cols <- grep("_SE$", names(df), value = TRUE)
  neff_cols <- grep("_NEFF$", names(df), value = TRUE)

  # Initialize the result data frame with the specified outcome columns
  result_df <- data.frame(
    beta.Y = df[[paste0(outcome_trait, "_BETA")]],
    se.Y = df[[paste0(outcome_trait, "_SE")]],
    neff.Y = df[[paste0(outcome_trait, "_NEFF")]]
  )
  names(result_df)[1:3] <- c(paste0("beta.", outcome_trait), paste0("se.", outcome_trait), paste0("N.", outcome_trait))

  # Initialize a list to store trait name mappings
  trait_mappings <- list()

  for (beta_col in beta_cols) {
    # Skip the outcome_BETA column
    if (beta_col == paste0(outcome_trait, "_BETA")) next

    # Get the corresponding SE and NEFF columns
    se_col <- gsub("_BETA$", "_SE", beta_col)
    neff_col <- gsub("_BETA$", "_NEFF", beta_col)

    # Extract the trait name from the column name
    trait_name <- gsub("_BETA$", "", beta_col)

    # Add the trait columns to the result data frame
    result_df[[paste0("beta.", trait_name)]] <- df[[beta_col]]
    result_df[[paste0("se.", trait_name)]] <- df[[se_col]]
    result_df[[paste0("N.", trait_name)]] <- df[[neff_col]]

    # Store the trait name mapping
    trait_mappings[[paste0("beta.", trait_name)]] <- beta_col
    trait_mappings[[paste0("se.", trait_name)]] <- se_col
    trait_mappings[[paste0("N.", trait_name)]] <- neff_col
  }

  # Standardize BETA and SE for all traits
  for (trait in names(trait_mappings)) {
    if (grepl("beta", trait)) {
      trait_name <- gsub("beta.", "", trait)
      c <- result_df[[paste0("se.", trait_name)]] * sqrt(result_df[[paste0("N.", trait_name)]])
      result_df[[paste0("beta.", trait_name)]] <- result_df[[paste0("beta.", trait_name)]] / c
      result_df[[paste0("se.", trait_name)]] <- result_df[[paste0("se.", trait_name)]] / c
    }
  }

  result_df <- na.omit(result_df)

  return(result_df)
}


