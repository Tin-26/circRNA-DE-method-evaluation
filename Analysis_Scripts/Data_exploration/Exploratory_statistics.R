library(dplyr)
library(edgeR)
library(ggplot2)

datasets <- list(
  `CIRI_HCC-Tissue` = "/media/meteor/FatDawg/Benchmark_Paper/PRJNA716508_hcc_study2/circRNA_outs/CIRI3/HCC-tissue_CIRI-Candidate.BSJ_Matrix",
  `CE2_HCC-Tissue`   = "/media/meteor/FatDawg/Benchmark_Paper/PRJNA716508_hcc_study2/circRNA_outs/CE2/HCC-tissue_CE2-BSJ_Matrix.txt",
  `CIRI_BC-Tissue` = "/media/meteor/FatDawg/Benchmark_Paper/PRJNA553624_breast_cancer/circRNA_outs/CIRI3/BC_CIRI-Candidate.BSJ_Matrix",
  `CE2_BC-Tissue`  = "/media/meteor/FatDawg/Benchmark_Paper/PRJNA553624_breast_cancer/circRNA_outs/CE2/BC_CE2-BSJ_Matrix.txt",
  `CIRI_HCC-PBMC` = "/media/meteor/ChunkyBoi/BC-alarm/CircRNA/PBMC/38samples/circRNA_outs/CIRI3/HCC-PBMC_CIRI-Candidate.BSJ_Matrix",
  `CE2_HCC-PBMC`  = "/media/meteor/ChunkyBoi/BC-alarm/CircRNA/PBMC/38samples/circRNA_outs/CE2/HCC-PBMC_CE2-BSJ_Matrix.txt",
  `CIRI_EBC1` = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2026-02-19/circRNA_outs/CIRI3/CIRI-Candidate.BSJ_Matrix",
  `CE2_EBC1` = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2026-02-19/circRNA_outs/CE2/CE2-BSJ_Matrix.txt",
  `CIRI_EBC2` = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-25/circRNA_outs/CIRI3/EBC2_CIRI-Candidate.BSJ_Matrix",
  `CE2_EBC2` = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-25/circRNA_outs/CE2/EBC2_CE2-BSJ_Matrix.txt"
)

metadata_paths <- list(
  `HCC-Tissue` = "/media/meteor/FatDawg/Benchmark_Paper/metadata_online_datasets/PRJNA716508_metadata_cleaned.csv",
  `BC-Tissue` = "/media/meteor/FatDawg/Benchmark_Paper/metadata_online_datasets/PRJNA553624_metadata.csv",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/metadata_online_datasets/PRJNA754685_metadata_cleaned.csv",
  `EBC1` = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2026-02-19/Metadata.csv",
  `EBC2` = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-25/11-25-Metadata.csv"
  )

filter_labels <- c(
  "Unfiltered"  = "Unfiltered",
  "default"     = "AutoFilter",
  "Default"     = "AutoFilter",
  "5"           = "Min 5",
  "min.count=5" = "Min 5",
  "1"           = "Min 1",
  "min.count=1" = "Min 1",
  "Pre-filter"  = "Pre-filter"
)

results <- list()
raw_ids  <- list()
raw_matched <- list()
raw_matrices <- list()

for (name in names(datasets)) {
  
  df <- read.table(datasets[[name]], sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  annotation_col <- c("gene_id", "geneName")
  df <- df[, !(colnames(df) %in% annotation_col), drop = FALSE]
  
  if (grepl("HCC-Tissue", name)) meta <- read.csv(metadata_paths$`HCC-Tissue`)
  if (grepl("BC-Tissue",  name)) meta <- read.csv(metadata_paths$`BC-Tissue`)
  if (grepl("HCC-PBMC",   name)) meta <- read.csv(metadata_paths$`HCC-PBMC`)
  if (grepl("EBC1",       name)) meta <- read.csv(metadata_paths$`EBC1`)
  if (grepl("EBC2",       name)) meta <- read.csv(metadata_paths$`EBC2`)
  
  meta_sub <- meta %>% filter(Sample.Name %in% colnames(df))
  # Store raw IDs before any filtering
  raw_ids[[name]] <- rownames(df)
  raw_matrices[[name]] <- df
  # Compute pre-filter matched for CIRI-CE2 pairs
  # Only do this when we hit a CIRI entry and its CE2 partner is already loaded
  if (grepl("^CIRI_", name)) {
    ce2_name <- sub("^CIRI_", "CE2_", name)
    if (!is.null(raw_ids[[ce2_name]])) {
      mapped <- normalize_and_match_ids(raw_ids[[name]], raw_ids[[ce2_name]])
      raw_matched[[name]] <- sum(!is.na(mapped))
      cat(sprintf("Pre-filter matched for %s: %d\n", name, raw_matched[[name]]))
    }
  } else if (grepl("^CE2_", name)) {
    ciri_name <- sub("^CE2_", "CIRI_", name)
    if (!is.null(raw_ids[[ciri_name]])) {
      mapped <- normalize_and_match_ids(raw_ids[[ciri_name]], raw_ids[[name]])
      raw_matched[[ciri_name]] <- sum(!is.na(mapped))
      cat(sprintf("Pre-filter matched for %s: %d\n", ciri_name, raw_matched[[ciri_name]]))
    }
  }
  
  # --- Per-sample stats BEFORE filtering (computed once, shared across filter versions) ---
  median_before              <- apply(df, 2, median)
  average_before             <- apply(df, 2, mean)
  percent_zero_before        <- apply(df, 2, function(x) mean(x == 0) * 100)
  circRNAs_per_sample_before <- apply(df, 2, function(x) sum(x > 0))
  bsj0_before <- apply(df, 2, function(x) sum(x == 0))
  bsj1_before <- apply(df, 2, function(x) sum(x == 1))
  bsj2_before <- apply(df, 2, function(x) sum(x == 2))
  bsj3_before <- apply(df, 2, function(x) sum(x == 3))
  
  for (mc in c("default", "1", "5")) {
    if (mc == "default") {
      keep <- filterByExpr(df, group = meta_sub$Group)
    } else {
      keep <- filterByExpr(df, group = meta_sub$Group, min.count = as.numeric(mc))
    }
    filtered_df <- df[keep, , drop = FALSE]
    
    # --- Library sizes ---
    dge <- DGEList(counts = filtered_df, group = meta_sub$Group)
    dge <- calcNormFactors(dge, method = "TMM")
    
    eff_lib_size <- dge$samples$lib.size * dge$samples$norm.factors
    
    lib_size_df <- data.frame(
      Sample         = rownames(dge$samples),
      raw_lib_size   = dge$samples$lib.size,
      norm_factor    = dge$samples$norm.factors,
      eff_lib_size   = eff_lib_size,
      Group          = meta_sub$Group[match(rownames(dge$samples), meta_sub$Sample.Name)],
      row.names      = NULL
    )
    
    cat(name, "|", mc, "-> Retained", nrow(filtered_df), "circRNAs\n")
    
    # --- Per-sample stats AFTER filtering ---
    median_after              <- apply(filtered_df, 2, median)
    average_after             <- apply(filtered_df, 2, mean)
    percent_zero_after        <- apply(filtered_df, 2, function(x) mean(x == 0) * 100)
    circRNAs_per_sample_after <- apply(filtered_df, 2, function(x) sum(x > 0))
    bsj0_after <- apply(filtered_df, 2, function(x) sum(x == 0))
    bsj1_after <- apply(filtered_df, 2, function(x) sum(x == 1))
    bsj2_after <- apply(filtered_df, 2, function(x) sum(x == 2))
    bsj3_after <- apply(filtered_df, 2, function(x) sum(x == 3))
    
    summary_table <- data.frame(
      Sample                 = colnames(df),
      circRNAs_before        = circRNAs_per_sample_before,
      circRNAs_after         = circRNAs_per_sample_after,
      Median_before          = median_before,
      Median_after           = median_after,
      Average_before         = average_before,
      Average_after          = average_after,
      Percent_zero_before    = percent_zero_before,
      Percent_zero_after     = percent_zero_after,
      BSJ0_before            = bsj0_before,
      BSJ0_after             = bsj0_after,
      BSJ1_before            = bsj1_before,
      BSJ1_after             = bsj1_after,
      BSJ2_before            = bsj2_before,
      BSJ2_after             = bsj2_after,
      BSJ3_before            = bsj3_before,
      BSJ3_after             = bsj3_after,
      row.names = NULL,
      check.names = FALSE
    )
    
    # Correlations (unchanged logic)
    corr_by_cond <- list()
    for (cond in unique(meta_sub$Group)) {
      samples_cond <- intersect(meta_sub$Sample.Name[meta_sub$Group == cond], colnames(df))
      if (length(samples_cond) > 1) {
        mat_cond <- df[, samples_cond, drop = FALSE]
        n <- length(samples_cond)
        corr_s <- pval_s <- stat_s <- matrix(NA, n, n)
        corr_p <- pval_p <- stat_p <- matrix(NA, n, n)
        rownames(corr_s) <- colnames(corr_s) <- samples_cond
        rownames(pval_s) <- colnames(pval_s) <- samples_cond
        rownames(stat_s) <- colnames(stat_s) <- samples_cond
        rownames(corr_p) <- colnames(corr_p) <- samples_cond
        rownames(pval_p) <- colnames(pval_p) <- samples_cond
        rownames(stat_p) <- colnames(stat_p) <- samples_cond
        for (i in 1:n) {
          for (j in i:n) {
            x <- mat_cond[, i]; y <- mat_cond[, j]
            xj <- x + runif(length(x), 0, 1e-6)
            yj <- y + runif(length(y), 0, 1e-6)
            ts <- cor.test(xj, yj, method = "spearman")
            corr_s[i,j] <- corr_s[j,i] <- ts$estimate
            pval_s[i,j] <- pval_s[j,i] <- ts$p.value
            stat_s[i,j] <- stat_s[j,i] <- ts$statistic
            tp <- cor.test(x, y, method = "pearson")
            corr_p[i,j] <- corr_p[j,i] <- tp$estimate
            pval_p[i,j] <- pval_p[j,i] <- tp$p.value
            stat_p[i,j] <- stat_p[j,i] <- tp$statistic
          }
        }
        corr_by_cond[[cond]] <- list(
          spearman = list(correlation = corr_s, pvalue = pval_s, statistic = stat_s),
          pearson  = list(correlation = corr_p, pvalue = pval_p, statistic = stat_p)
        )
      } else {
        cat("Condition", cond, "has <2 samples, skipping.\n")
      }
    }
    
    # Store under a combined key e.g. "CIRI_HCC-Tissue_default"
    result_key <- paste(name, mc, sep = "_")
    results[[result_key]] <- list(
      summary_table   = summary_table,
      correlations    = corr_by_cond,
      filtered_matrix = filtered_df,
      lib_sizes = lib_size_df
    )
  }
}

output_dir <- "/media/meteor/FatDawg/Benchmark_Paper/Dataset_statistics"

write_if_new <- function(x, path, row.names = TRUE) {
  if (file.exists(path)) {
    cat("  Skipping (exists):", basename(path), "\n")
    return(invisible(NULL))
  }
  write.csv(x, path, row.names = row.names)
}

for (dataset_name in names(results)) {
  
  res <- results[[dataset_name]]
  
  # Extract filter version from key and create subdirectory
  filter_version <- str_extract(dataset_name, "default|1|5")
  sub_dir        <- file.path(output_dir, filter_version)
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## ---------- SUMMARY ----------
  write_if_new(
    res$summary_table,
    file.path(sub_dir, paste0(dataset_name, "_summary_table.csv")),
    row.names = FALSE
  )
  
  ## ---------- CORRELATIONS ----------
  for (cond in names(res$correlations)) {
    for (method in c("spearman", "pearson")) {
      
      corr_list <- res$correlations[[cond]][[method]]
      if (is.null(corr_list)) next
      
      write_if_new(corr_list$correlation,
                   file.path(sub_dir, paste0(dataset_name, "_", cond, "_", method, "_correlation.csv")))
      write_if_new(corr_list$pvalue,
                   file.path(sub_dir, paste0(dataset_name, "_", cond, "_", method, "_pvalue.csv")))
      write_if_new(corr_list$statistic,
                   file.path(sub_dir, paste0(dataset_name, "_", cond, "_", method, "_statistic.csv")))
    }
  }
  
  cat("Done:", dataset_name, "\n")
}

long_summary <- function(summary_df, dataset_name, identifier, filter_version) {
  summary_df %>%
    mutate(
      Dataset        = dataset_name,
      Identifier     = identifier,
      Filter_version = filter_version,
      Sample         = Sample
    ) %>%
    tidyr::pivot_longer(
      cols      = -c(Sample, Dataset, Identifier, Filter_version),
      names_to  = "Metric",
      values_to = "Value"
    )
}

# Combine multiple datasets
library(stringr)

summary_stats <- dplyr::bind_rows(lapply(names(results), function(key) {
  
  parts      <- str_match(key, "^(CIRI|CE2)_(.+)_(default|1|5)$")
  identifier <- parts[, 2]   # CIRI or CE2
  dataset    <- parts[, 3]   # HCC-Tissue, BC-Tissue etc
  filter_ver <- parts[, 4]   # default, 1, 5
  
  tbl <- results[[key]]$summary_table
  if (is.null(tbl)) return(NULL)
  
  data.frame(
    identifier              = identifier,
    dataset                 = dataset,
    filter_version          = filter_ver,
    avg_pct_zero_before     = round(mean(tbl$Percent_zero_before, na.rm = TRUE), 2),
    avg_pct_zero_after      = round(mean(tbl$Percent_zero_after,  na.rm = TRUE), 2),
    avg_bsj_count_after     = round(mean(tbl$Average_after,       na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
}))

# Pivot wide so CIRI and CE2 are side by side
summary_wide <- summary_stats %>%
  tidyr::pivot_wider(
    id_cols     = c(dataset, filter_version),
    names_from  = identifier,
    values_from = c(avg_pct_zero_before, avg_pct_zero_after, avg_bsj_count_after)
  ) %>%
  dplyr::select(
    Dataset                         = dataset,
    Filter                          = filter_version,
    `CIRI3 Zero counts pre-filter (%)` = avg_pct_zero_before_CIRI,
    `CE2 Zero counts pre-filter (%)`   = avg_pct_zero_before_CE2,
    `CIRI3 post-filter avg BSJ`        = avg_bsj_count_after_CIRI,
    `CE2 post-filter avg BSJ`          = avg_bsj_count_after_CE2,
    `CIRI3 post-filter zero (%)`       = avg_pct_zero_after_CIRI,
    `CE2 post-filter zero (%)`         = avg_pct_zero_after_CE2
  ) %>%
  dplyr::mutate(Filter = factor(Filter, levels = c("default", "1", "5"))) %>%
  dplyr::arrange(Filter, Dataset)

print(summary_wide)
