#### Pipeline to run different DE tools and normalization methods on multiple datasets ####

# Authors: Erda Qorri & Valentin Varga

####################
# Import libraries #
####################
library(edgeR)
library(DESeq2)
library(limma)

library(tidyverse)
library(tibble)

########
# Seed #
########
set.seed(42)

###### DE Functions ######

#### Simulated data ####

##################
### Limma-voom ###
##################

## ZeroSet ##
# Functions #
run_limma_default_DE0 <- function(bsj_counts_path,
                                  metadata_path,
                                  contrast_formula = NULL,
                                  edgeR_norm = "TMM", # switch to TMMwsp and quantile
                                  pvalue_thresholds = c(0.01, 0.05, 0.10),
                                  rep_id = NA) {
  # Start timer
  set.seed(42)
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    
    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    bsj_fc_md$sample_id = colnames(bsj_fc)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
    }
    
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    group <- factor(bsj_fc_md$Group)
    
    print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    cat(sprintf("✓ Detected groups: %s\n", paste(levels(group), collapse = " vs ")))
    cat(sprintf("✓ Sample sizes: %s\n", paste(table(group), collapse = " and ")))
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 3: Normalization #
    cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
    
    dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
    
    # Step 4: Design matrix and model fitting
    cat("Starting Step 4: Fitting voomLmFit model... \n")
    
    # Design matrix (no intercept was included since it is a simple two group comparison)
    design_mtx = model.matrix(~ 0 + group)
    
    # Assigns the name of the levels to the design mtx colnames
    colnames(design_mtx) <- levels(group)
    
    # Sanity check
    
    if (is.null(contrast_formula)) {
      stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
    }
    
    # the provided colnames must reference the existing design columns
    cn <- colnames(design_mtx)
    
    parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
    
    if (!all(parts %in% cn)) {
      stop(sprintf(
        "Contrast groups not found in design. Provided: %s. Available: %s",
        paste(parts, collapse = ", "),
        paste(cn, collapse = ", ")
      ))
    }
    
    contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
    print(contrasts)
    
    #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
    # it replaces all the commands with one function
    dgelist_v <- voom(dgelist_norm, design_mtx, plot = F)
    
    dgelist_vfit <- lmFit(dgelist_v, design_mtx)
    
    # Step 5: Set the contrast matrix and DE analysis
    cat(" Starting Step 5: Testing contrasts... \n")
    
    # Get the group names from design matrix columns
    cat(sprintf(" Testing contrast: %s\n", contrast_formula))
    
    dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
    
    dgelist_efit <- eBayes(dgelist_vfit)
    
    # Extract the full results
    summary(decideTests(dgelist_efit))
    de_results <- topTable(dgelist_efit, number = Inf, sort.by="none")
    
    # Step 6: Calculate the metrics at each Pvalue threshold #
    cat(" Starting Step 6: Calculating metrics at each Pvalue threshold... \n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    
    #Sig genes
    sig_genes_list <- list()
    # For each FDR threshold
    for (pval in pvalue_thresholds) {
      n_sig = sum(de_results$P.Value < pval, na.rm = TRUE)
      #n_sig = sum(results$PValue < pval, na.rm = TRUE)
      # false positive rate: number of genes that were found to be significant over the number of genes following the filtering step (genes tested for DE) 
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      print(pval_names)
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
      
      #Collect genes that were significant
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      
      # Get genes significant at this threshold
      sig_idx <- de_results$P.Value < pval
      sig_genes <- rownames(de_results)[sig_idx]
      
      # Store as comma-separated string
      sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(de_results$P.Value, na.rm = TRUE), 4)
    metrics$median_pval = round(median(de_results$P.Value, na.rm = TRUE), 4)
    metrics$min_pval = round(min(de_results$P.Value, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    
    return(list(metrics = metrics, de_results = de_results, sig_genes = sig_genes_list))
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(error_row)
  })
}

run_limma_DFRBT_DE0 <- function(bsj_counts_path,
                                metadata_path,
                                contrast_formula = NULL,
                                edgeR_norm = "TMM", # switch to TMMwsp
                                pvalue_thresholds = c(0.01, 0.05, 0.10),
                                rep_id = NA) {
  set.seed(42)
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    # Not an mtx -> converts it to one by removing the X col which is the sample id col
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    
    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    
    # Step 3: Add the sample id information to the bsj mtx for tracking purposes
    bsj_fc_md$sample_id = colnames(bsj_fc)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
    }
    
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    group <- factor(bsj_fc_md$Group)
    
    print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    cat(sprintf("Hey there! I detected the groups: %s\n", paste(levels(group), collapse = " vs ")))
    cat(sprintf("And the number of samples per group is: %s\n", paste(table(group), collapse = " and ")))
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    
    # pct of genes that passed the automatic filtering
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 3: Normalization #
    cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
    
    dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
    
    # Step 4: Design matrix and model fitting
    cat("Starting Step 4: Fitting voomLmFit model... \n")
    
    # Design matrix (no intercept was included since it is a simple two group comparison)
    design_mtx = model.matrix(~ 0 + group)
    
    # Assigns the name of the levels to the design mtx colnames
    colnames(design_mtx) <- levels(group)
    
    # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
    
    if (is.null(contrast_formula)) {
      stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
    }
    
    # the provided colnames must reference the existing design columns
    cn <- colnames(design_mtx)
    
    parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
    
    if (!all(parts %in% cn)) {
      stop(sprintf(
        "Contrast groups not found in design. Provided: %s. Available: %s",
        paste(parts, collapse = ", "),
        paste(cn, collapse = ", ")
      ))
    }
    
    contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
    print(contrasts)
    
    #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
    # it replaces all the commands with one function
    
    # Fit the default limma pipeline but using quantile normalization
    dgelist_v <- voom(dgelist_norm, design_mtx)
    
    # Fit the original linear model with lmFit just set robust to TRUE
    dgelist_vfit <- lmFit(dgelist_v, design_mtx, robust = TRUE)
    
    # Step 5: Set the contrast matrix and DE analysis
    cat(" Starting Step 5: Testing contrasts... \n")
    
    # Get the group names from design matrix columns
    cat(sprintf(" Testing contrast: %s\n", contrast_formula))
    
    dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
    
    dgelist_efit <- eBayes(dgelist_vfit)
    
    # Extract the full results
    summary(decideTests(dgelist_efit))
    de_results <- topTable(dgelist_efit, number = Inf, sort.by="none")
    
    # Step 6: Calculate the metrics at each Pvalue threshold #
    cat(" Starting Step 6: Calculating metrics at each Pvalue threshold... \n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    
    sig_genes_list <- list()
    # For each FDR threshold
    for (pval in pvalue_thresholds) {
      n_sig = sum(de_results$P.Value < pval, na.rm = TRUE)
      #n_sig = sum(results$PValue < pval, na.rm = TRUE)
      # false positive rate: number of genes that were found to be significant over the number of genes following the filtering step (genes tested for DE) 
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      print(pval_names)
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
      
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      
      # Get genes significant at this threshold
      sig_idx <- de_results$P.Value < pval
      sig_genes <- rownames(de_results)[sig_idx]
      
      # Store as comma-separated string
      sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(de_results$P.Value, na.rm = TRUE), 4)
    metrics$median_pval = round(median(de_results$P.Value, na.rm = TRUE), 4)
    metrics$min_pval = round(min(de_results$P.Value, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    
    return(list(metrics = metrics, de_results = de_results, sig_genes = sig_genes_list ))
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(error_row)
  })
}

run_limma_voomlmFit_DE0 <- function(bsj_counts_path,
                                    metadata_path,
                                    contrast_formula = NULL,
                                    edgeR_norm = "TMM", # switch to TMMwsp and quantile
                                    pvalue_thresholds = c(0.01, 0.05, 0.10),
                                    rep_id = NA) {
  # Start timer
  set.seed(42)
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    
    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    bsj_fc_md$sample_id = colnames(bsj_fc)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    stopifnot("sample_id" %in% colnames(bsj_fc_md))
    
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), ]
    
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in the BSJ counts matrix and the metadata file do not match 1:1.")
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    group <- factor(bsj_fc_md$Group)
    
    print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    cat(sprintf("✓ Detected groups: %s\n", paste(levels(group), collapse = " vs ")))
    cat(sprintf("✓ Sample sizes: %s\n", paste(table(group), collapse = " and ")))
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 3: Normalization #
    cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
    
    dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
    
    # Step 4: Design matrix and model fitting
    cat("Starting Step 4: Fitting voomLmFit model... \n")
    
    # Design matrix (no intercept was included since it is a simple two group comparison)
    design_mtx = model.matrix(~ 0 + group)
    
    # Assigns the name of the levels to the design mtx colnames
    colnames(design_mtx) <- levels(group)
    
    # Sanity check
    
    if (is.null(contrast_formula)) {
      stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
    }
    
    # the provided colnames must reference the existing design columns
    cn <- colnames(design_mtx)
    
    parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
    
    if (!all(parts %in% cn)) {
      stop(sprintf(
        "Contrast groups not found in design. Provided: %s. Available: %s",
        paste(parts, collapse = ", "),
        paste(cn, collapse = ", ")
      ))
    }
    
    contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
    print(contrasts)
    
    #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
    # it replaces all the commands with one function
    dgelist_v_fit <- voomLmFit(dgelist_norm, design = design_mtx, sample.weights = TRUE)
    
    # Step 5: Set the contrast matrix and DE analysis
    cat(" Starting Step 5: Testing contrasts... \n")
    
    # Get the group names from design matrix columns
    cat(sprintf(" Comparing: %s vs %s\n", cn[1], cn[2]))
    
    dgelist_v_fit <- contrasts.fit(dgelist_v_fit, contrasts = contrasts)
    
    dgelist_v_efit <- eBayes(dgelist_v_fit)
    
    # Extract the full results
    summary(decideTests(dgelist_v_efit))
    de_results <- topTable(dgelist_v_efit, number = Inf, sort.by="none")
    
    # Step 6: Calculate the metrics at each Pvalue threshold #
    cat(" Starting Step 6: Calculating metrics at each Pvalue threshold... \n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    
    sig_genes_list <- list()
    # For each FDR threshold
    for (pval in pvalue_thresholds) {
      
      # calculate the sum of the rows that fulfill the pval condition
      n_sig = sum(de_results$P.Value < pval, na.rm = TRUE)
      
      # false positive rate: number of genes that were found to be significant over the number of genes following the filtering step (genes tested for DE) 
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      print(pval_names)
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
      
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      
      # Get genes significant at this threshold
      sig_idx <- de_results$P.Value < pval
      sig_genes <- rownames(de_results)[sig_idx]
      
      # Store as comma-separated string
      sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(de_results$P.Value, na.rm = TRUE), 4)
    metrics$median_pval = round(median(de_results$P.Value, na.rm = TRUE), 4)
    metrics$min_pval = round(min(de_results$P.Value, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    
    # added a new option to retrieve the de_results as well in case are needed for validation
    return(list(metrics = metrics, de_results = de_results, sig_genes = sig_genes_list))
    # OUTSIDE the loop
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(error_row)
  })
}

run_limma_quantile_DE0 <- function(bsj_counts_path,
                                   metadata_path,
                                   contrast_formula = NULL,
                                   #edgeR_norm = none, # switch to TMMwsp
                                   pvalue_thresholds = c(0.01, 0.05, 0.10),
                                   rep_id = NA) {
  set.seed(42)
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    # Not an mtx -> converts it to one by removing the X col which is the sample id col
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    
    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    
    # Step 3: Add the sample id information to the bsj mtx for tracking purposes
    bsj_fc_md$sample_id = colnames(bsj_fc)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
    }
    
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    group <- factor(bsj_fc_md$Group)
    
    print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    cat(sprintf("Hey there! I detected the groups: %s\n", paste(levels(group), collapse = " vs ")))
    cat(sprintf("And the number of samples per group is: %s\n", paste(table(group), collapse = " and ")))
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    
    # pct of genes that passed the automatic filtering
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 3: Normalization #
    #cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
    
    #dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
    
    # Step 4: Design matrix and model fitting
    cat("Starting Step 4: Fitting voomLmFit model... \n")
    
    # Design matrix (no intercept was included since it is a simple two group comparison)
    design_mtx = model.matrix(~ 0 + group)
    
    # Assigns the name of the levels to the design mtx colnames
    colnames(design_mtx) <- levels(group)
    
    # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
    
    if (is.null(contrast_formula)) {
      stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
    }
    
    # the provided colnames must reference the existing design columns
    cn <- colnames(design_mtx)
    
    parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
    
    if (!all(parts %in% cn)) {
      stop(sprintf(
        "Contrast groups not found in design. Provided: %s. Available: %s",
        paste(parts, collapse = ", "),
        paste(cn, collapse = ", ")
      ))
    }
    
    contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
    print(contrasts)
    
    #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
    # it replaces all the commands with one function
    
    # Fit the default limma pipeline but using quantile normalization
    dgelist_v <- voom(dgelist_filtered, design_mtx, normalize.method = "quantile")
    
    # Fit the original linear model with lmFit just set robust to TRUE
    dgelist_vfit <- lmFit(dgelist_v, design_mtx)
    
    # Step 5: Set the contrast matrix and DE analysis
    cat(" Starting Step 5: Testing contrasts... \n")
    
    # Get the group names from design matrix columns
    cat(sprintf(" Testing contrast: %s\n", contrast_formula))
    
    dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
    
    dgelist_efit <- eBayes(dgelist_vfit)
    
    # Extract the full results
    summary(decideTests(dgelist_efit))
    de_results <- topTable(dgelist_efit, number = Inf, sort.by="none")
    
    # Step 6: Calculate the metrics at each Pvalue threshold #
    cat(" Starting Step 6: Calculating metrics at each Pvalue threshold... \n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    
    sig_genes_list <- list()
    # For each FDR threshold
    for (pval in pvalue_thresholds) {
      n_sig = sum(de_results$P.Value < pval, na.rm = TRUE)
      #n_sig = sum(results$PValue < pval, na.rm = TRUE)
      # false positive rate: number of genes that were found to be significant over the number of genes following the filtering step (genes tested for DE) 
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      print(pval_names)
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
      
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      
      # Get genes significant at this threshold
      sig_idx <- de_results$P.Value < pval
      sig_genes <- rownames(de_results)[sig_idx]
      # Store as comma-separated string
      sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(de_results$P.Value, na.rm = TRUE), 4)
    metrics$median_pval = round(median(de_results$P.Value, na.rm = TRUE), 4)
    metrics$min_pval = round(min(de_results$P.Value, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    
    return(list(metrics = metrics, de_results = de_results, sig_genes = sig_genes_list))
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(error_row)
  })
}

# Iterators #
limma_voom_quantile <- function(data_dir,
                                dataset_name,
                                file_pattern,
                                #norm_method = "TMM",
                                contrast_formula,
                                pvalue_thresholds = c(0.01, 0.05, 0.10),
                                output_dir = NULL,
                                filter = NULL,
                                de_label=NULL,
                                variant= NULL) {
  set.seed(42)
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  #cat(sprintf("Normalization method: %s\n", norm_method))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  tp_files <- gsub("_counts\\.csv", "_TP.csv", all_files)
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    tp_path <- paste0(data_dir, tp_files[i])
    
    # Run analysis
    if (de_label=="DE0"){
      result <- run_limma_quantile_DE0(
        bsj_counts_path = counts_path,
        metadata_path = design_path,
        # edgeR_norm = norm_method,
        contrast_formula = contrast_formula,
        pvalue_thresholds = pvalue_thresholds,
        rep_id = rep_numbers[i]
      )
      sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), result$sig_genes)
      
    } else {
      result <- run_limma_vquantile_with_truth(
        bsj_counts_path = counts_path,
        metadata_path = design_path,
        truth_set_path = tp_path,
        # edgeR_norm = norm_method,
        contrast_formula = contrast_formula,
        pvalue_thresholds = pvalue_thresholds,
        rep_id = rep_numbers[i])
    }
    # Store result
    results_list[[i]] <- result$metrics
    
    # Track success/failure
    
  }
  
  #### Step 5: Combine results ####
  cat("Combining results...\n")
  results_df <- do.call(rbind, results_list)
  
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      norm_method = "Quantile",
      .before = 1
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_%s_%s_Limma-voom_quantile_summary.csv", dataset_name, filter, de_label)
  output_path <- file.path(output_dir, output_filename)
  
  write.csv(results_df, output_path, row.names = FALSE)
  
  sig_genes_filename <- sprintf("%s_%s_%s_Limma-voom_quantile_sig_genes.csv", 
                                dataset_name, filter, de_label)
  sig_genes_path <- file.path(output_dir, sig_genes_filename)
  write.csv(sig_genes_df, sig_genes_path, row.names = FALSE)
  
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(results_df)
}

limma_voom_DFRBT <- function(data_dir,
                             dataset_name,
                             file_pattern,
                             norm_method = "TMM",
                             contrast_formula,
                             pvalue_thresholds = c(0.01, 0.05, 0.10),
                             output_dir = NULL,
                             filter = NULL,
                             de_label = NULL,
                             variant= NULL) {
  set.seed(42)
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat(sprintf("Normalization method: %s\n", norm_method))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    # Run analysis
    if (de_label == "DE0") { 
      result <- run_limma_DFRBT_DE0(
        bsj_counts_path = counts_path,
        metadata_path = design_path,
        edgeR_norm = norm_method,
        contrast_formula = contrast_formula,
        pvalue_thresholds = pvalue_thresholds,
        rep_id = rep_numbers[i]
        )
      sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), result$sig_genes)
      
    } else {
      result <- run_limma_vDFBRT_with_truth(
        bsj_counts_path = counts_path,
        metadata_path = design_path,
        edgeR_norm = norm_method,
        contrast_formula = contrast_formula,
        pvalue_thresholds = pvalue_thresholds,
        rep_id = rep_numbers[i]
      )
    }
    # Store result
    results_list[[i]] <- result$metrics
    
    # Track success/failure
    
  }
  
  #### Step 5: Combine results ####
  cat("Combining results...\n")
  results_df <- do.call(rbind, results_list)
  
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      norm_method = norm_method,
      .before = 1
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_%s_DE0_Limma-voom_%s_%s_summary.csv", dataset_name, filter, variant, norm_method)
  output_path <- file.path(output_dir, output_filename)
  
  write.csv(results_df, output_path, row.names = FALSE)
  sig_genes_filename <- sprintf("%s_%s_DE0_Limma-voom_%s_%s_sig_genes.csv", 
                                dataset_name, filter, variant, norm_method)
  write.csv(sig_genes_df, file.path(output_dir, sig_genes_filename), row.names = FALSE)
  
  
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(list(metrics = results_df, sig_genes = sig_genes_df))
}

limma_voom_TMM_default <- function(data_dir,
                                   dataset_name,
                                   file_pattern,
                                   norm_method = "TMM",
                                   contrast_formula,
                                   pvalue_thresholds = c(0.01, 0.05, 0.10),
                                   output_dir = NULL,
                                   filter=NULL,
                                   de_label = NULL,
                                   variant= NULL) {
  set.seed(42)
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat(sprintf("Normalization method: %s\n", norm_method))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    # Run analysis
    if (de_label == "DE0") {
      result <- run_limma_default_DE0(
        bsj_counts_path = counts_path,
        metadata_path = design_path,
        edgeR_norm = norm_method,
        contrast_formula = contrast_formula,
        pvalue_thresholds = pvalue_thresholds,
        rep_id = rep_numbers[i]
      )
      sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), result$sig_genes)
    
    } else {
      result <- run_limma_vdefault_with_truth(
        bsj_counts_path = counts_path,
        metadata_path = design_path,
        edgeR_norm = norm_method,
        contrast_formula = contrast_formula,
        pvalue_thresholds = pvalue_thresholds,
        rep_id = rep_numbers[i])
    }
    # Store result
    results_list[[i]] <- result$metrics
    
    # Track success/failure
    
  }
  
  #### Step 5: Combine results ####
  cat("Combining results...\n")
  results_df <- do.call(rbind, results_list)
  
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      norm_method = norm_method,
      .before = 1
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_%s_DE0_Limma-voom_%s_%s_summary.csv", dataset_name, filter, variant, norm_method)
  output_path <- file.path(output_dir, output_filename)
  
  write.csv(results_df, output_path, row.names = FALSE)
  sig_genes_filename <- sprintf("%s_%s_DE0_Limma-voom_%s_%s_sig_genes.csv", 
                                dataset_name, filter, variant, norm_method)
  write.csv(sig_genes_df, file.path(output_dir, sig_genes_filename), row.names = FALSE)
  
  return(list(metrics = results_df, sig_genes = sig_genes_df))
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(results_df)
}

limma_voom_TMM_LmFit_weights <- function(data_dir,
                                         dataset_name,
                                         file_pattern,
                                         norm_method = "TMMwsp",
                                         contrast_formula = NULL,
                                         pvalue_thresholds = c(0.01, 0.05, 0.10),
                                         output_dir = NULL,
                                         filter = NULL,
                                         de_label = NULL,
                                         variant = NULL) {
  set.seed(42)
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat(sprintf("Normalization method: %s\n", norm_method))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) output_dir <- data_dir
  if (!endsWith(data_dir, "/")) data_dir <- paste0(data_dir, "/")
  
  #### Step 1: Discover count files ####
  cat("Discovering simulation files...\n")
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  if (length(all_files) == 0) stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and design files ####
  rep_numbers <- as.numeric(gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files))
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  if (!all(file.exists(paste0(data_dir, design_files)))) {
    missing <- design_files[!file.exists(paste0(data_dir, design_files))]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  #### Step 4: Process each replicate ####
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    result <- tryCatch({
      if (de_label == "DE0") {
        res <- run_limma_voomlmFit_DE0(
          bsj_counts_path = counts_path,
          metadata_path = design_path,
          edgeR_norm = norm_method,
          pvalue_thresholds = pvalue_thresholds,
          rep_id = rep_numbers[i],
          contrast_formula = contrast_formula
        )
        
        sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), res$sig_genes)
      } else {
        res <- run_limma_LmFit_with_truth(
          bsj_counts_path = counts_path,
          metadata_path = design_path,
          edgeR_norm = norm_method,
          pvalue_thresholds = pvalue_thresholds,
          rep_id = rep_numbers[i],
          contrast_formula = contrast_formula
        )
        # always safe
        sig_genes_list[[i]] <- data.frame(
          rep_id = rep_numbers[i],
          gene = character(0),
          stringsAsFactors = FALSE
        )
      }
      res$metrics
    }, error = function(e) {
      cat(sprintf("  ERROR in replicate %d: %s\n", rep_numbers[i], e$message))
      n_failed <<- n_failed + 1
      sig_genes_list[[i]] <<- data.frame(
        rep_id = rep_numbers[i],
        gene = character(0),
        stringsAsFactors = FALSE
      )
      NULL
    })
    
    results_list[[i]] <- result
    if (!is.null(result)) n_success <- n_success + 1
    cat("\n")
  }
  
  #### Step 5: Combine results ####
  cat("Combining results...\n")
  results_df <- do.call(rbind, results_list)
  
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  
  results_df <- results_df %>%
    mutate(dataset = dataset_name,
           norm_method = norm_method,
           .before = 1)
  
  #### Step 6: Save results ####
  output_filename <- sprintf("%s_%s_DE0_Limma-voom_%s_%s_summary.csv", dataset_name, filter, variant, norm_method)
  write.csv(results_df, file.path(output_dir, output_filename), row.names = FALSE)
  
  sig_genes_filename <- sprintf("%s_%s_DE0_Limma-voom_%s_%s_sig_genes.csv", dataset_name, filter, variant, norm_method)
  write.csv(sig_genes_df, file.path(output_dir, sig_genes_filename), row.names = FALSE)
  
  #### Step 7: Summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\nTotal replicates: %d\nSuccessful: %d\nFailed: %d\n", dataset_name, n_total, n_success, n_total - n_success))
  cat(sprintf("Results saved to:\n%s\n", file.path(output_dir, output_filename)))
  
  return(list(metrics = results_df, sig_genes = sig_genes_df))
}

## Signal Set ##
# Functions #
run_limma_vdefault_with_truth <- function(bsj_counts_path,
                                          metadata_path,
                                          truth_set_path, # path to the TP file from SPSimSeq
                                          contrast_formula, # unchanged
                                          edgeR_norm = "TMMwsp", # optional -> change to TMMwsp
                                          thresholds = c(0.01, 0.05, 0.10),
                                          threshold_col = c("adj.P.Val"),
                                          rep_id = NA) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  # Step 1: Load the data and the metadata #
  cat("Starting Step 1: Loading the BSJ matrix ...\n")
  
  bsj_fc <- read.csv(bsj_counts_path)
  rownames(bsj_fc) <- bsj_fc$X
  bsj_fc$X <- NULL
  
  # circRNAs before filtering with filterByExpr
  n_genes_initial_counts <- nrow(bsj_fc)
  
  # Load the metadata 
  bsj_fc_md <- read.csv(metadata_path)
  bsj_fc_md$sample_id <- colnames(bsj_fc)
  
  # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
  
  if (!"sample_id" %in% colnames(bsj_fc_md)) {
    stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
  }
  
  bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
  
  if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
    stop("Sample IDs in counts and metadata do not match 1:1.")
  }
  
  ##### Truth set loading #####
  cat("Starting Step 2: Loading the Truth Set ...\n")
  
  truth_set <- read.csv(truth_set_path)
  
  # Check that the colnames X and DE.ind are present in the truth set
  if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
    stop("Truth file must contain columns: X, DE.ind")
  }
  
  # to ensure this is logical (True and False values)
  truth_set$DE.ind <- as.logical(truth_set$DE.ind)
  
  # set the number of circRNAs that were in the truth set
  n_genes_initial_truth_set <- length(unique(truth_set$X))
  
  ##### DGEList and filtering #####
  cat("Starting Step 3: Creating the DGEList and starting the automated filtering ...\n")
  
  bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
  group <- factor(bsj_fc_md$Group)
  
  # must always be only.2 groups for comparison in our case
  if (length(levels(group)) != 2) stop("Expected exactly 2 groups.")
  
  dgelist <- DGEList(counts = bsj_fc, group = group)
  keep <- filterByExpr(dgelist, group = group)
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # extract the universe size -> circRNA ids that will be tested by the DE methods
  tested_circRNAs <- rownames(dgelist_filtered)
  n_circRNAs_tested <- length(tested_circRNAs)
  
  # calculate the pct of retained circRNAs post filtering
  pct_retained <- round(n_circRNAs_tested / n_genes_initial_truth_set * 100, 2)
  
  # Restrict truth to TESTED universe
  truth_set_tested <- truth_set %>% filter(X %in% tested_circRNAs)
  
  # If truth has any duplicates, remove them to avoid join inflation
  truth_set_tested <- truth_set_tested %>% distinct(X, .keep_all = TRUE)
  
  ##### Normalizing and setting the contrasts #####
  cat("Starting Step 4: Normalizing the DGEList and setting the contrasts as defined by Alexa ...\n")  
  
  dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = edgeR_norm)
  
  # Design matrix (no intercept was included since it is a simple two group comparison)
  design_mtx = model.matrix(~ 0 + group)
  
  # Assigns the name of the levels to the design mtx colnames
  colnames(design_mtx) <- levels(group)
  
  # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
  
  if (is.null(contrast_formula)) {
    stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
  }
  
  # the provided colnames must reference the existing design columns
  cn <- colnames(design_mtx)
  
  parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
  
  if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
  
  if (!all(parts %in% cn)) {
    stop(sprintf(
      "Contrast groups not found in design. Provided: %s. Available: %s",
      paste(parts, collapse = ", "),
      paste(cn, collapse = ", ")
    ))
  }
  
  contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
  print(contrasts)
  
  ##### Fit glmQLFit model and test with glmQLFTest the contrasts #####
  cat("Starting Step 5: Fitting voom and lmFit and testing the contrasts ...\n")  
  #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
  # it replaces all the commands with one function
  dgelist_v <- voom(dgelist_filtered_norm, design_mtx, plot = F)
  
  dgelist_vfit <- lmFit(dgelist_v, design_mtx)
  
  # Step 5: Set the contrast matrix and DE analysis
  cat(" Starting Step 5: Testing contrasts... \n")
  
  # Get the group names from design matrix columns
  cat(sprintf(" Testing contrast: %s\n", contrast_formula))
  
  dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
  
  dgelist_efit <- eBayes(dgelist_vfit)
  
  # Extract the full results
  summary(decideTests(dgelist_efit))
  de_results <- topTable(dgelist_efit, number = Inf, sort.by="none")
  
  ##### Create the results table #####
  cat("Starting Step 6: Calculating the metrics and generating the results table ...\n")  
  de_results <- as.data.frame(de_results)
  de_results$X <- rownames(de_results)
  rownames(de_results) <- NULL
  
  # Ensure we only evaluate on tested IDs that appear in results
  # (should match, but we keep this safe)
  tested_ids_try2 <- de_results$X
  truth_tested2 <- truth_set_tested %>% dplyr::filter(X %in% tested_ids_try2)
  truth_tested2 <- truth_tested2 %>% dplyr::distinct(X, .keep_all = TRUE)
  
  # Extract the true DE (TP) and non DE (TN)
  DE_truth <- truth_tested2 %>% dplyr::filter(DE.ind == TRUE)
  DE_notTruth <- truth_tested2 %>% dplyr::filter(DE.ind == FALSE)
  
  # Helper for safe division
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  
  # Calculate metrics for every specified threshold
  out <- lapply(thresholds, function(th) {
    
    called <- de_results %>% filter(.data[[threshold_col]] <= th)
    not_called <- de_results %>% filter(.data[[threshold_col]] > th)
    
    TP <- nrow(inner_join(called,     DE_truth,  by = "X"))
    FP <- nrow(inner_join(called,     DE_notTruth, by = "X"))
    TN <- nrow(inner_join(not_called, DE_notTruth, by = "X"))
    FN <- nrow(inner_join(not_called, DE_truth,  by = "X"))
    
    precision <- safe_div(TP, TP + FP)
    recall    <- safe_div(TP, TP + FN)  # = sensitivity = TPR
    specificity <- safe_div(TN, TN + FP)
    fpr <- safe_div(FP, FP + TN)
    fdr <- safe_div(FP, TP + FP)
    f1  <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                  NA_real_,
                  2 * precision * recall / (precision + recall))
    
    data.frame(
      rep_id = rep_id,
      threshold_type = threshold_col,
      threshold = th,
      
      initial_genes_truth = n_genes_initial_truth_set,
      initial_genes_counts = n_genes_initial_counts,
      tested_genes = n_circRNAs_tested,
      pct_retained = pct_retained,
      
      TP = TP, FP = FP, TN = TN, FN = FN,
      
      sensitivity = recall * 100,
      specificity = specificity * 100,
      precision = precision,
      recall = recall,
      F1 = f1,
      FDR = fdr * 100,
      FPR = fpr * 100,
      
      n_called = nrow(called)
    )
  }) %>% bind_rows()
  
  runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  out$runtime_sec <- runtime_sec
  
  list(metrics = out, de_results = de_results)
}

run_limma_vDFBRT_with_truth <- function(bsj_counts_path,
                                        metadata_path,
                                        truth_set_path, # path to the TP file from SPSimSeq
                                        contrast_formula, # unchanged
                                        edgeR_norm = "TMM", # optional -> change to TMM
                                        thresholds = c(0.01, 0.05, 0.10),
                                        threshold_col = c("adj.P.Val"),
                                        rep_id = NA) {
  
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  # Step 1: Load the data and the metadata #
  cat("Starting Step 1: Loading the BSJ matrix ...\n")
  
  bsj_fc <- read.csv(bsj_counts_path)
  rownames(bsj_fc) <- bsj_fc$X
  bsj_fc$X <- NULL
  
  # circRNAs before filtering with filterByExpr
  n_genes_initial_counts <- nrow(bsj_fc)
  
  # Load the metadata 
  bsj_fc_md <- read.csv(metadata_path)
  bsj_fc_md$sample_id <- colnames(bsj_fc)
  
  # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
  
  if (!"sample_id" %in% colnames(bsj_fc_md)) {
    stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
  }
  
  bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
  
  if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
    stop("Sample IDs in counts and metadata do not match 1:1.")
  }
  
  ##### Truth set loading #####
  cat("Starting Step 2: Loading the Truth Set ...\n")
  
  truth_set <- read.csv(truth_set_path)
  
  # Check that the colnames X and DE.ind are present in the truth set
  if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
    stop("Truth file must contain columns: X, DE.ind")
  }
  
  # to ensure this is logical (True and False values)
  truth_set$DE.ind <- as.logical(truth_set$DE.ind)
  
  # set the number of circRNAs that were in the truth set
  n_genes_initial_truth_set <- length(unique(truth_set$X))
  
  ##### DGEList and filtering #####
  cat("Starting Step 3: Creating the DGEList and starting the automated filtering ...\n")
  
  bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
  group <- factor(bsj_fc_md$Group)
  
  # must always be only.2 groups for comparison in our case
  if (length(levels(group)) != 2) stop("Expected exactly 2 groups.")
  
  dgelist <- DGEList(counts = bsj_fc, group = group)
  keep <- filterByExpr(dgelist, group = group)
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # extract the universe size -> circRNA ids that will be tested by the DE methods
  tested_circRNAs <- rownames(dgelist_filtered)
  n_circRNAs_tested <- length(tested_circRNAs)
  
  # calculate the pct of retained circRNAs post filtering
  pct_retained <- round(n_circRNAs_tested / n_genes_initial_truth_set * 100, 2)
  
  # Restrict truth to TESTED universe
  truth_set_tested <- truth_set %>% filter(X %in% tested_circRNAs)
  
  # If truth has any duplicates, remove them to avoid join inflation
  truth_set_tested <- truth_set_tested %>% distinct(X, .keep_all = TRUE)
  
  ##### Normalizing and setting the contrasts #####
  cat("Starting Step 4: Normalizing the DGEList and setting the contrasts as defined by Alexa ...\n")  
  
  dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = edgeR_norm)
  
  # Design matrix (no intercept was included since it is a simple two group comparison)
  design_mtx = model.matrix(~ 0 + group)
  
  # Assigns the name of the levels to the design mtx colnames
  colnames(design_mtx) <- levels(group)
  
  # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
  
  if (is.null(contrast_formula)) {
    stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
  }
  
  # the provided colnames must reference the existing design columns
  cn <- colnames(design_mtx)
  
  parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
  
  if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
  
  if (!all(parts %in% cn)) {
    stop(sprintf(
      "Contrast groups not found in design. Provided: %s. Available: %s",
      paste(parts, collapse = ", "),
      paste(cn, collapse = ", ")
    ))
  }
  
  contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
  print(contrasts)
  
  ##### Fit glmQLFit model and test with glmQLFTest the contrasts #####
  cat("Starting Step 5: Fitting voom and lmFit and testing the contrasts ...\n")  
  #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
  # it replaces all the commands with one function
  dgelist_v <- voom(dgelist_filtered_norm, design_mtx)
  
  dgelist_vfit <- lmFit(dgelist_v, design_mtx, robust = TRUE)
  
  # Step 5: Set the contrast matrix and DE analysis
  cat(" Starting Step 5: Testing contrasts... \n")
  
  # Get the group names from design matrix columns
  cat(sprintf(" Testing contrast: %s\n", contrast_formula))
  
  dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
  
  dgelist_efit <- eBayes(dgelist_vfit)
  
  # Extract the full results
  summary(decideTests(dgelist_efit))
  de_results <- topTable(dgelist_efit, number = Inf)
  
  ##### Create the results table #####
  cat("Starting Step 6: Calculating the metrics and generating the results table ...\n")  
  de_results <- as.data.frame(de_results)
  de_results$X <- rownames(de_results)
  rownames(de_results) <- NULL
  
  # Ensure we only evaluate on tested IDs that appear in results
  # (should match, but we keep this safe)
  tested_ids_try2 <- de_results$X
  truth_tested2 <- truth_set_tested %>% dplyr::filter(X %in% tested_ids_try2)
  truth_tested2 <- truth_tested2 %>% dplyr::distinct(X, .keep_all = TRUE)
  
  # Extract the true DE (TP) and non DE (TN)
  DE_truth <- truth_tested2 %>% dplyr::filter(DE.ind == TRUE)
  DE_notTruth <- truth_tested2 %>% dplyr::filter(DE.ind == FALSE)
  
  # Helper for safe division
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  
  # Calculate metrics for every specified threshold
  out <- lapply(thresholds, function(th) {
    
    called <- de_results %>% filter(.data[[threshold_col]] <= th)
    not_called <- de_results %>% filter(!(X %in% called$X))
    
    TP <- nrow(inner_join(called,     DE_truth,  by = "X"))
    FP <- nrow(inner_join(called,     DE_notTruth, by = "X"))
    TN <- nrow(inner_join(not_called, DE_notTruth, by = "X"))
    FN <- nrow(inner_join(not_called, DE_truth,  by = "X"))
    
    precision <- safe_div(TP, TP + FP)
    recall    <- safe_div(TP, TP + FN)  # = sensitivity = TPR
    specificity <- safe_div(TN, TN + FP)
    fpr <- safe_div(FP, FP + TN)
    fdr <- safe_div(FP, TP + FP)
    f1  <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                  NA_real_,
                  2 * precision * recall / (precision + recall))
    
    data.frame(
      rep_id = rep_id,
      threshold_type = threshold_col,
      threshold = th,
      
      initial_genes_truth = n_genes_initial_truth_set,
      initial_genes_counts = n_genes_initial_counts,
      tested_genes = n_circRNAs_tested,
      pct_retained = pct_retained,
      
      TP = TP, FP = FP, TN = TN, FN = FN,
      
      sensitivity = recall * 100,
      specificity = specificity * 100,
      precision = precision,
      recall = recall,
      F1 = f1,
      FDR = fdr * 100,
      FPR = fpr * 100,
      
      n_called = nrow(called)
    )
  }) %>% bind_rows()
  
  runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  out$runtime_sec <- runtime_sec
  
  list(metrics = out, de_results = de_results)
}

run_limma_LmFit_with_truth <- function(bsj_counts_path,
                                       metadata_path,
                                       truth_set_path, # path to the TP file from SPSimSeq
                                       contrast_formula, # unchanged
                                       edgeR_norm = "TMMwsp", # optional -> change to TMMwsp
                                       thresholds = c(0.01, 0.05, 0.10),
                                       threshold_col = c("adj.P.Val"),
                                       rep_id = NA) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  # Step 1: Load the data and the metadata #
  cat("Starting Step 1: Loading the BSJ matrix ...\n")
  
  bsj_fc <- read.csv(bsj_counts_path)
  rownames(bsj_fc) <- bsj_fc$X
  bsj_fc$X <- NULL
  
  # circRNAs before filtering with filterByExpr
  n_genes_initial_counts <- nrow(bsj_fc)
  
  # Load the metadata 
  bsj_fc_md <- read.csv(metadata_path)
  bsj_fc_md$sample_id <- colnames(bsj_fc)
  
  # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
  
  if (!"sample_id" %in% colnames(bsj_fc_md)) {
    stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
  }
  
  bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
  
  if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
    stop("Sample IDs in counts and metadata do not match 1:1.")
  }
  
  ##### Truth set loading #####
  cat("Starting Step 2: Loading the Truth Set ...\n")
  
  truth_set <- read.csv(truth_set_path)
  
  # Check that the colnames X and DE.ind are present in the truth set
  if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
    stop("Truth file must contain columns: X, DE.ind")
  }
  
  # to ensure this is logical (True and False values)
  truth_set$DE.ind <- as.logical(truth_set$DE.ind)
  
  # set the number of circRNAs that were in the truth set
  n_genes_initial_truth_set <- length(unique(truth_set$X))
  
  ##### DGEList and filtering #####
  cat("Starting Step 3: Creating the DGEList and starting the automated filtering ...\n")
  
  bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
  group <- factor(bsj_fc_md$Group)
  
  # must always be only.2 groups for comparison in our case
  if (length(levels(group)) != 2) stop("Expected exactly 2 groups.")
  
  dgelist <- DGEList(counts = bsj_fc, group = group)
  keep <- filterByExpr(dgelist, group = group)
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # extract the universe size -> circRNA ids that will be tested by the DE methods
  tested_circRNAs <- rownames(dgelist_filtered)
  n_circRNAs_tested <- length(tested_circRNAs)
  
  # calculate the pct of retained circRNAs post filtering
  pct_retained <- round(n_circRNAs_tested / n_genes_initial_truth_set * 100, 2)
  
  # Restrict truth to TESTED universe
  truth_set_tested <- truth_set %>% filter(X %in% tested_circRNAs)
  
  # If truth has any duplicates, remove them to avoid join inflation
  truth_set_tested <- truth_set_tested %>% distinct(X, .keep_all = TRUE)
  
  ##### Normalizing and setting the contrasts #####
  cat("Starting Step 4: Normalizing the DGEList and setting the contrasts as defined by Alexa ...\n")  
  
  dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = edgeR_norm)
  
  # Design matrix (no intercept was included since it is a simple two group comparison)
  design_mtx = model.matrix(~ 0 + group)
  
  # Assigns the name of the levels to the design mtx colnames
  colnames(design_mtx) <- levels(group)
  
  # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
  
  if (is.null(contrast_formula)) {
    stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
  }
  
  # the provided colnames must reference the existing design columns
  cn <- colnames(design_mtx)
  
  parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
  
  if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
  
  if (!all(parts %in% cn)) {
    stop(sprintf(
      "Contrast groups not found in design. Provided: %s. Available: %s",
      paste(parts, collapse = ", "),
      paste(cn, collapse = ", ")
    ))
  }
  
  contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
  print(contrasts)
  
  ##### Fit glmQLFit model and test with glmQLFTest the contrasts #####
  cat("Starting Step 5: Fitting voom and lmFit and testing the contrasts ...\n")  
  #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
  # it replaces all the commands with one function
  dgelist_v_fit <- voomLmFit(dgelist_filtered_norm, design = design_mtx, sample.weights = TRUE)
  
  # Step 5: Set the contrast matrix and DE analysis
  cat(" Starting Step 5: Testing contrasts... \n")
  
  # Get the group names from design matrix columns
  cat(sprintf(" Testing contrast: %s\n", contrast_formula))
  
  dgelist_vfit <- contrasts.fit(dgelist_v_fit, contrasts = contrasts)
  
  dgelist_efit <- eBayes(dgelist_vfit)
  
  # Extract the full results
  summary(decideTests(dgelist_efit))
  de_results <- topTable(dgelist_efit, number = Inf)
  
  ##### Create the results table #####
  cat("Starting Step 6: Calculating the metrics and generating the results table ...\n")  
  de_results <- as.data.frame(de_results)
  de_results$X <- rownames(de_results)
  rownames(de_results) <- NULL
  
  # Ensure we only evaluate on tested IDs that appear in results
  # (should match, but we keep this safe)
  tested_ids_try2 <- de_results$X
  truth_tested2 <- truth_set_tested %>% dplyr::filter(X %in% tested_ids_try2)
  truth_tested2 <- truth_tested2 %>% dplyr::distinct(X, .keep_all = TRUE)
  
  # Extract the true DE (TP) and non DE (TN)
  DE_truth <- truth_tested2 %>% dplyr::filter(DE.ind == TRUE)
  DE_notTruth <- truth_tested2 %>% dplyr::filter(DE.ind == FALSE)
  
  # Helper for safe division
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  
  # Calculate metrics for every specified threshold
  out <- lapply(thresholds, function(th) {
    
    called <- de_results %>% filter(.data[[threshold_col]] <= th)
    not_called <- de_results %>% filter(!(X %in% called$X))
    
    TP <- nrow(inner_join(called,     DE_truth,  by = "X"))
    FP <- nrow(inner_join(called,     DE_notTruth, by = "X"))
    TN <- nrow(inner_join(not_called, DE_notTruth, by = "X"))
    FN <- nrow(inner_join(not_called, DE_truth,  by = "X"))
    
    precision <- safe_div(TP, TP + FP)
    recall    <- safe_div(TP, TP + FN)  # = sensitivity = TPR
    specificity <- safe_div(TN, TN + FP)
    fpr <- safe_div(FP, FP + TN)
    fdr <- safe_div(FP, TP + FP)
    f1  <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                  NA_real_,
                  2 * precision * recall / (precision + recall))
    
    data.frame(
      rep_id = rep_id,
      threshold_type = threshold_col,
      threshold = th,
      
      initial_genes_truth = n_genes_initial_truth_set,
      initial_genes_counts = n_genes_initial_counts,
      tested_genes = n_circRNAs_tested,
      pct_retained = pct_retained,
      
      TP = TP, FP = FP, TN = TN, FN = FN,
      
      sensitivity = recall * 100,
      specificity = specificity * 100,
      precision = precision,
      recall = recall,
      F1 = f1,
      FDR = fdr * 100,
      FPR = fpr * 100,
      
      n_called = nrow(called)
    )
  }) %>% bind_rows()
  
  runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  out$runtime_sec <- runtime_sec
  
  list(metrics = out, de_results = de_results)
}

run_limma_vquantile_with_truth <- function(bsj_counts_path,
                                           metadata_path,
                                           truth_set_path, # path to the TP file from SPSimSeq
                                           contrast_formula, # unchanged
                                           # edgeR_norm = NULL, 
                                           thresholds = c(0.01, 0.05, 0.10),
                                           threshold_col = "adj.P.Val",
                                           rep_id = NA) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  # Step 1: Load the data and the metadata #
  cat("Starting Step 1: Loading the BSJ matrix ...\n")
  
  bsj_fc <- read.csv(bsj_counts_path)
  rownames(bsj_fc) <- bsj_fc$X
  bsj_fc$X <- NULL
  
  # circRNAs before filtering with filterByExpr
  n_genes_initial_counts <- nrow(bsj_fc)
  
  # Load the metadata 
  bsj_fc_md <- read.csv(metadata_path)
  bsj_fc_md$sample_id <- colnames(bsj_fc)
  
  # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
  
  if (!"sample_id" %in% colnames(bsj_fc_md)) {
    stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
  }
  
  bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
  
  if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
    stop("Sample IDs in counts and metadata do not match 1:1.")
  }
  
  ##### Truth set loading #####
  cat("Starting Step 2: Loading the Truth Set ...\n")
  
  truth_set <- read.csv(truth_set_path)
  
  # Check that the colnames X and DE.ind are present in the truth set
  if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
    stop("Truth file must contain columns: X, DE.ind")
  }
  
  # to ensure this is logical (True and False values)
  truth_set$DE.ind <- as.logical(truth_set$DE.ind)
  
  # set the number of circRNAs that were in the truth set
  n_genes_initial_truth_set <- length(unique(truth_set$X))
  
  ##### DGEList and filtering #####
  cat("Starting Step 3: Creating the DGEList and starting the automated filtering ...\n")
  
  bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
  group <- factor(bsj_fc_md$Group)
  
  # must always be only.2 groups for comparison in our case
  if (length(levels(group)) != 2) stop("Expected exactly 2 groups.")
  
  dgelist <- DGEList(counts = bsj_fc, group = group)
  keep <- filterByExpr(dgelist, group = group)
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # extract the universe size -> circRNA ids that will be tested by the DE methods
  tested_circRNAs <- rownames(dgelist_filtered)
  n_circRNAs_tested <- length(tested_circRNAs)
  
  # calculate the pct of retained circRNAs post filtering
  pct_retained <- round(n_circRNAs_tested / n_genes_initial_truth_set * 100, 2)
  
  # Restrict truth to TESTED universe
  truth_set_tested <- truth_set %>% filter(X %in% tested_circRNAs)
  
  # If truth has any duplicates, remove them to avoid join inflation
  truth_set_tested <- truth_set_tested %>% distinct(X, .keep_all = TRUE)
  
  ##### Normalizing and setting the contrasts #####
  cat("Starting Step 4: Normalizing the DGEList and setting the contrasts as defined by Alexa ...\n")  
  
  #dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = edgeR_norm)
  
  # Design matrix (no intercept was included since it is a simple two group comparison)
  design_mtx = model.matrix(~ 0 + group)
  
  # Assigns the name of the levels to the design mtx colnames
  colnames(design_mtx) <- levels(group)
  
  # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
  
  if (is.null(contrast_formula)) {
    stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
  }
  
  # the provided colnames must reference the existing design columns
  cn <- colnames(design_mtx)
  
  parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
  
  if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
  
  if (!all(parts %in% cn)) {
    stop(sprintf(
      "Contrast groups not found in design. Provided: %s. Available: %s",
      paste(parts, collapse = ", "),
      paste(cn, collapse = ", ")
    ))
  }
  
  contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
  print(contrasts)
  
  ##### Fit glmQLFit model and test with glmQLFTest the contrasts #####
  cat("Starting Step 5: Fitting voom and lmFit and testing the contrasts ...\n")
  
  #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
  # it replaces all the commands with one function
  dgelist_v <- voom(dgelist_filtered, design_mtx, normalize.method = "quantile")
  
  dgelist_vfit <- lmFit(dgelist_v, design_mtx)
  
  # Step 5: Set the contrast matrix and DE analysis
  cat(" Starting Step 5: Testing contrasts... \n")
  
  # Get the group names from design matrix columns
  cat(sprintf(" Testing contrast: %s\n", contrast_formula))
  
  dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
  
  dgelist_efit <- eBayes(dgelist_vfit)
  
  # Extract the full results
  summary(decideTests(dgelist_efit))
  de_results <- topTable(dgelist_efit, number = Inf)
  
  ##### Create the results table #####
  cat("Starting Step 6: Calculating the metrics and generating the results table ...\n")  
  de_results <- as.data.frame(de_results)
  de_results$X <- rownames(de_results)
  rownames(de_results) <- NULL
  
  # Ensure we only evaluate on tested IDs that appear in results
  # (should match, but we keep this safe)
  tested_ids_try2 <- de_results$X
  truth_tested2 <- truth_set_tested %>% dplyr::filter(X %in% tested_ids_try2)
  truth_tested2 <- truth_tested2 %>% dplyr::distinct(X, .keep_all = TRUE)
  
  # Extract the true DE (TP) and non DE (TN)
  DE_truth <- truth_tested2 %>% dplyr::filter(DE.ind == TRUE)
  DE_notTruth <- truth_tested2 %>% dplyr::filter(DE.ind == FALSE)
  
  # Helper for safe division
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  
  # Calculate metrics for every specified threshold
  out <- lapply(thresholds, function(th) {
    
    called <- de_results %>% filter(.data[[threshold_col]] <= th)
    not_called <- de_results %>% filter(!(X %in% called$X))
    
    TP <- nrow(inner_join(called,     DE_truth,  by = "X"))
    FP <- nrow(inner_join(called,     DE_notTruth, by = "X"))
    TN <- nrow(inner_join(not_called, DE_notTruth, by = "X"))
    FN <- nrow(inner_join(not_called, DE_truth,  by = "X"))
    
    precision <- safe_div(TP, TP + FP)
    recall    <- safe_div(TP, TP + FN)  # = sensitivity = TPR
    specificity <- safe_div(TN, TN + FP)
    fpr <- safe_div(FP, FP + TN)
    fdr <- safe_div(FP, TP + FP)
    f1  <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                  NA_real_,
                  2 * precision * recall / (precision + recall))
    
    data.frame(
      rep_id = rep_id,
      threshold_type = threshold_col,
      threshold = th,
      
      initial_genes_truth = n_genes_initial_truth_set,
      initial_genes_counts = n_genes_initial_counts,
      tested_genes = n_circRNAs_tested,
      pct_retained = pct_retained,
      
      TP = TP, FP = FP, TN = TN, FN = FN,
      
      sensitivity = recall * 100,
      specificity = specificity * 100,
      precision = precision,
      recall = recall,
      F1 = f1,
      FDR = fdr * 100,
      FPR = fpr * 100,
      
      n_called = nrow(called)
    )
  }) %>% bind_rows()
  
  runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  out$runtime_sec <- runtime_sec
  
  list(metrics = out, de_results = de_results)
}

# Special iterator #
lv_quantile_signal <- function(data_dir,
                               dataset_name,
                               file_pattern,
                               runner_fun,               # <- pass run_limma_voom_with_truth or run_edgeR_with_truth
                               #  norm_method = "TMM",
                               contrast_formula,
                               thresholds = c(0.01, 0.05, 0.10),
                               threshold_col = "adj.P.Val",
                               output_dir = NULL,
                               filter_type=NULL) {
  set.seed(42)  
  threshold_col <- match.arg(threshold_col)
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  #cat(sprintf("Normalization method: %s\n", norm_method))
  cat(sprintf("Threshold column: %s\n", threshold_col))
  cat("========================================\n\n")
  
  if (is.null(output_dir)) output_dir <- data_dir
  if (!endsWith(data_dir, "/")) data_dir <- paste0(data_dir, "/")
  
  ##### Step 1: Find all count files #####
  cat("Discovering simulation files...\n")
  count_files <- list.files(
    data_dir,
    pattern = paste0(file_pattern, ".*_counts\\.csv$"),
    full.names = FALSE
  )
  
  if (length(count_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  cat(sprintf("Found %d count files\n\n", length(count_files)))
  
  ##### Step 2: Extract dataset rep numbers and the count files #####
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", count_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Order the number of the replicates
  ord <- order(rep_numbers)
  cat("Ordering replicates:\n")
  cat(sprintf("  Original order: %s\n", paste(rep_numbers, collapse = ", ")))
  
  count_files <- count_files[ord]
  rep_numbers <- rep_numbers[ord]
  cat(sprintf("  Sorted order: %s\n\n", paste(rep_numbers, collapse = ", ")))
  
  # Extract the design and truth files
  design_files <- gsub("_counts\\.csv$", "_design.csv", count_files)
  truth_files  <- gsub("_counts\\.csv$", "_TP.csv", count_files)
  
  # Check existence of these files in the defined data dir path
  design_exist <- file.exists(paste0(data_dir, design_files))
  truth_exist  <- file.exists(paste0(data_dir, truth_files))
  
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  if (!all(truth_exist)) {
    missing <- truth_files[!truth_exist]
    stop(sprintf("Missing TP truth files: %s", paste(missing, collapse = ", ")))
  }
  
  ##### Step 3: Initiate result storing in a list #####
  results_list <- vector("list", length(count_files))
  de_results_list <- vector("list", length(count_files))
  # as long as the number of count files
  n_total <- length(count_files)
  n_success <- 0
  n_failed <- 0
  
  #### Step 4: Process each replicate dataset #####
  for (i in seq_along(count_files)) {
    rep_id <- rep_numbers[i] # gets updated with each iteration
    cat("\n----------------------------------------\n")
    cat(sprintf("Processing replicate %d/%d (rep_%d)\n", i, n_total, rep_id))
    cat("----------------------------------------\n")
    
    counts_path <- paste0(data_dir, count_files[i])
    design_path <- paste0(data_dir, design_files[i])
    truth_path  <- paste0(data_dir, truth_files[i])
    
    # Sanity Check 1: Print which files are being loaded for this iteration
    cat(sprintf("Files for rep_%d:\n", rep_id))
    cat(sprintf("  Counts:  %s\n", count_files[i]))
    cat(sprintf("  Design:  %s\n", design_files[i]))
    cat(sprintf("  Truth:   %s\n", truth_files[i]))
    
    # Sanity Check 2: Confirm files exist and check file sizes
    cat("\nFile verification:\n")
    cat(sprintf("  Counts exists: %s (size: %d bytes)\n", 
                file.exists(counts_path), 
                ifelse(file.exists(counts_path), file.info(counts_path)$size, 0)))
    cat(sprintf("  Design exists: %s (size: %d bytes)\n", 
                file.exists(design_path), 
                ifelse(file.exists(design_path), file.info(design_path)$size, 0)))
    cat(sprintf("  Truth exists:  %s (size: %d bytes)\n", 
                file.exists(truth_path), 
                ifelse(file.exists(truth_path), file.info(truth_path)$size, 0)))
    
    # Sanity Check 3: Extract rep number from each filename to ensure they match
    counts_rep <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", count_files[i])
    design_rep <- gsub(paste0(".*", file_pattern, "(\\d+)_design\\.csv"), "\\1", design_files[i])
    truth_rep  <- gsub(paste0(".*", file_pattern, "(\\d+)_TP\\.csv"), "\\1", truth_files[i])
    
    cat(sprintf("\nReplicate ID consistency check:\n"))
    cat(sprintf("  Expected rep_id: %d\n", rep_id))
    cat(sprintf("  Counts file rep: %s\n", counts_rep))
    cat(sprintf("  Design file rep: %s\n", design_rep))
    cat(sprintf("  Truth file rep:  %s\n", truth_rep))
    
    if (counts_rep != design_rep || counts_rep != truth_rep) {
      warning(sprintf("MISMATCH: Files for iteration %d do not all correspond to the same replicate!", i))
    } else {
      cat(sprintf(" ATTENTION: All files correspond to rep_%d\n", rep_id))
    }
    
    cat("\nRunning analysis...\n")
    
    # write the results file
    res <- tryCatch({
      runner_fun(
        bsj_counts_path  = counts_path,
        metadata_path    = design_path,
        truth_set_path   = truth_path,
        contrast_formula = contrast_formula,
        #  edgeR_norm       = norm_method,   # runner may ignore if not edgeR; keep for consistency
        thresholds       = thresholds,
        threshold_col    = threshold_col,
        rep_id           = rep_id
      )
    }, error = function(e) {
      cat(sprintf(" WARNING: ERROR in rep_%d: %s\n", rep_id, e$message))
      return(list(metrics = data.frame(
        rep_id = rep_id,
        error_message = e$message,
        stringsAsFactors = FALSE
      )))
    })
    
    # Store metrics (expected: data.frame with one row per threshold)
    results_list[[i]] <- res$metrics
    de_results_list[[i]] <- res$de_results
    # Update counters
    if ("error_message" %in% colnames(res$metrics)) {
      n_failed <- n_failed + 1
      cat(sprintf("  Status: FAILED\n"))
    } else {
      n_success <- n_success + 1
      cat(sprintf("  Status: SUCCESS\n"))
    }
  }
  
  ##### Step 5: Combine the results per threshold level ####
  cat("\n========================================\n")
  cat("Combining results...\n")
  results_df <- dplyr::bind_rows(results_list)
  names(de_results_list) <- paste0("rep_", rep_numbers)
  ##### Step 6: Add metadata  #####
  results_df <- results_df %>%
    dplyr::mutate(
      dataset = dataset_name,
      norm_method = "quantile",
      .before = 1
    )
  
  ##### Step 7: Save the output file #####
  output_filename <- sprintf("%s_%s_DE0.1_%s_%s_summary.csv",
                             dataset_name,
                             filter_type,
                             "Limma-voom",
                             "DFRBT_quantile")
  output_path <- file.path(output_dir, output_filename)
  write.csv(results_df, output_path, row.names = FALSE)
  
  de_results_combined <- dplyr::bind_rows(
    lapply(names(de_results_list), function(rep_name) {
      if (!is.null(de_results_list[[rep_name]])) {
        de_results_list[[rep_name]]$rep_id <- rep_name
        return(de_results_list[[rep_name]])
      }
    })
  )
  de_output_filename <- gsub("_summary.csv", "_DE_results_all_reps.csv", output_filename)
  de_output_path <- file.path(output_dir, de_output_filename)
  write.csv(de_results_combined, de_output_path, row.names = FALSE)
  
  
  ##### Step 8: Print summary #####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  cat(sprintf("\nDE_Results saved to:\n%s\n", de_output_path))
  # Some quick summaries (only if those columns exist)
  if (n_success > 0 && "tested_genes" %in% colnames(results_df)) {
    # Filter out failed replicates if error_message column exists
    if ("error_message" %in% colnames(results_df)) {
      valid <- results_df %>% dplyr::filter(is.na(error_message))
    } else {
      valid <- results_df
    }
    cat("\n--- Summary Statistics (tested universe) ---\n")
    cat(sprintf("Mean genes tested after filtering: %.0f (±%.0f)\n",
                mean(valid$tested_genes, na.rm = TRUE),
                sd(valid$tested_genes, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  return(list(
    metrics    = results_df,
    de_results = de_results_list
  ))
}

#############
### edgeR ###
#############

## ZeroSet ##
# Function #
run_edgeR_DE0 <- function(bsj_counts_path,
                          metadata_path,
                          edgeR_norm = "TMM",
                          pvalue_thresholds = c(0.01, 0.05, 0.10),
                          rep_id = NA,
                          contrast_formula = NULL) {
  set.seed(42)
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    #head(bsj_fc)

    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    bsj_fc_md$sample_id = colnames(bsj_fc)
    #head(bsj_fc_md)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    if(!identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      # Try to reorder the metadata to match count columns
      bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), ]
      
      # Verify again
      if(!identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
        stop("Sample IDs in the BSJ counts matrix and the metadata files do not match!")
      }
    
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    group = as.factor(bsj_fc_md$Group)
    #print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    #head(dgelist)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    #head(dgelist_filtered)
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    #print(n_genes_tested)
    
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    #print(n_genes_initial)
    #print(pct_filtered)
    
  cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
              n_genes_tested, n_genes_initial, pct_filtered))
  
  # Step 3: Normalization #
  cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
  
  dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
  #head(dgelist_norm)
  
  # Step 4: Design matrix and model fitting
  cat("Starting Step 4: Fitting glmQLFit model... \n")
  
  # Create valid R names for the design matrix groups
  group = factor(group)
  levels(group) = make.names(levels(group))
  #print(group)
  
  # Design matrix (no intercept was included since it is a simple two group comparison)
  design_mtx = model.matrix(~ 0 + group)
  colnames(design_mtx) = levels(group)
  #print(design_mtx)
  
  # Fit the quasi-likelihood model
  fit = glmQLFit(dgelist_norm, design = design_mtx, robust = TRUE)
  #dim(dgelist_norm) %>% print()
  
  # Step 5: Set the contrast matrix and DE analysis
  cat(" Starting Step 5: Testing contrasts... \n")
  
  # Create contrast: cancer vs healthy controls
  # contrast_formula = paste(levels(group)[1], "-", levels(group)[2])
  #print(contrast_formula)
  contrasts = makeContrasts(contrasts = contrast_formula, levels = design_mtx)
  #print(contrasts)
  
  # Quasi-likelihood F-test
  qlf = glmQLFTest(fit, contrast = contrasts)
  #print(qlf)
  
  # Extract the full results
  results = as.data.frame(topTags(qlf, n = Inf))
  # write.csv(results, "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/res.csv")

  # Step 6: Calculate the metrics at each Pvalue threshold #
  cat(" Starting Step 6: Calculating metrics at each Pvalue threshold... \n")
  
  metrics = data.frame(rep_id = rep_id)
  metrics$n_genes_initial = n_genes_initial
  metrics$n_genes_tested = n_genes_tested
  metrics$pct_filtered = pct_filtered
  sig_genes_list <- list()
  
  # For each FDR threshold
  for (pval in pvalue_thresholds) {
    n_sig = sum(results$PValue < pval, na.rm = TRUE)
    #n_sig = sum(results$PValue < pval, na.rm = TRUE)
    # false positive rate: number of genes that were found to be significant over the number of genes following the filtering step (genes tested for DE) 
    fpr = n_sig / n_genes_tested

    # Create column names based on the threshold
    pval_names = gsub("\\.", "", sprintf("%.3f", pval))
    print(pval_names)
    
    metrics[[paste0("n_sig_", pval_names)]] <- n_sig
    metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
    
    pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
    
    # Get genes significant at this threshold
    sig_idx <- results$PValue < pval
    sig_genes <- rownames(results)[sig_idx]
    
    # Store as comma-separated string
    sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    }
  
  # P-value distribution metrics (OUTSIDE the loop)
  metrics$mean_pval = round(mean(results$PValue, na.rm = TRUE), 4)
  metrics$median_pval = round(median(results$PValue, na.rm = TRUE), 4)
  metrics$min_pval = round(min(results$PValue, na.rm = TRUE), 6)
  
  # Runtime (OUTSIDE the loop)
  processing_time <- Sys.time()
  metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
  
  cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))
  
  return(list(
    metrics = metrics, 
    de_results = results,
    sig_genes = sig_genes_list  
  )) # OUTSIDE the loop
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(list(
      metrics = error_metrics,
      de_results = NULL,
      sig_genes = list()
    ))
  })
}

# Iterator #
analyze_edger_DE0 <- function(data_dir,
                                          dataset_name,
                                          file_pattern,
                                          norm_method = "TMM",
                                          pvalue_thresholds = c(0.01, 0.05, 0.10),
                                          output_dir = NULL,
                                          contrast_formula = NULL,
                              filter=NULL) {
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat(sprintf("Normalization method: %s\n", norm_method))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    # Run analysis
    result <- run_edgeR_DE0(
      bsj_counts_path = counts_path,
      metadata_path = design_path,
      edgeR_norm = norm_method,
      pvalue_thresholds = pvalue_thresholds,
      rep_id = rep_numbers[i],
      contrast_formula = contrast_formula
    )
    
    sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), result$sig_genes)
    # Store result
    results_list[[i]] <- result$metrics
    
    
    # Track success/failure
    if (is.na(result$metrics$n_genes_tested)) {
      n_failed <- n_failed + 1
    } else {
      n_success <- n_success + 1
    }
    
    cat("\n")
  }
  
  #### Step 5: Combine results ####
  cat("Combining results...\n")
  results_df <- dplyr::bind_rows(results_list)
  results_df <- as.data.frame(results_df)
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      norm_method = norm_method,
      .before = 1
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_%s_DE0_edgeR_%s_summary.csv", dataset_name, filter, norm_method)
  output_path <- file.path(output_dir, output_filename)
  sig_genes_filename <- sprintf("%s_%s_DE0_edgeRm_%s_sig_genes.csv", dataset_name, filter, norm_method)
  write.csv(sig_genes_df, file.path(output_dir, sig_genes_filename), row.names = FALSE)
  write.csv(results_df, output_path, row.names = FALSE)
  
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(results_df)
}

## Signal Set ##
# Function #
run_edgeR_vdefault_with_truth <- function(bsj_counts_path,
                                          metadata_path,
                                          truth_set_path, # path to the TP file from SPSimSeq
                                          contrast_formula, # unchanged
                                          edgeR_norm = "TMM", 
                                          thresholds = c(0.01, 0.05, 0.10),
                                          threshold_col = "FDR",
                                          rep_id = NA) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  # Step 1: Load the data and the metadata #
  cat("Starting Step 1: Loading the BSJ matrix ...\n")
  
  bsj_fc <- read.csv(bsj_counts_path)
  rownames(bsj_fc) <- bsj_fc$X
  bsj_fc$X <- NULL
  
  # circRNAs before filtering with filterByExpr
  n_genes_initial_counts <- nrow(bsj_fc)
  
  # Load the metadata 
  bsj_fc_md <- read.csv(metadata_path)
  bsj_fc_md$sample_id <- colnames(bsj_fc)
  
  # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
  
  if (!"sample_id" %in% colnames(bsj_fc_md)) {
    stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
  }
  
  bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
  
  if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
    stop("Sample IDs in counts and metadata do not match 1:1.")
  }
  
  ##### Truth set loading #####
  cat("Starting Step 2: Loading the Truth Set ...\n")
  
  truth_set <- read.csv(truth_set_path)
  
  # Check that the colnames X and DE.ind are present in the truth set
  if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
    stop("Truth file must contain columns: X, DE.ind")
  }
  
  # to ensure this is logical (True and False values)
  truth_set$DE.ind <- as.logical(truth_set$DE.ind)
  
  # set the number of circRNAs that were in the truth set
  n_genes_initial_truth_set <- length(unique(truth_set$X))
  
  ##### DGEList and filtering #####
  cat("Starting Step 3: Creating the DGEList and starting the automated filtering ...\n")
  
  bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
  group <- factor(bsj_fc_md$Group)
  
  # must always be only.2 groups for comparison in our case
  if (length(levels(group)) != 2) stop("Expected exactly 2 groups.")
  
  dgelist <- DGEList(counts = bsj_fc, group = group)
  keep <- filterByExpr(dgelist, group = group)
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # extract the universe size -> circRNA ids that will be tested by the DE methods
  tested_circRNAs <- rownames(dgelist_filtered)
  n_circRNAs_tested <- length(tested_circRNAs)
  
  # calculate the pct of retained circRNAs post filtering
  pct_retained <- round(n_circRNAs_tested / n_genes_initial_truth_set * 100, 2)
  
  # Restrict truth to TESTED universe
  truth_set_tested <- truth_set %>% filter(X %in% tested_circRNAs)
  
  # If truth has any duplicates, remove them to avoid join inflation
  truth_set_tested <- truth_set_tested %>% distinct(X, .keep_all = TRUE)
  
  ##### Normalizing and setting the contrasts #####
  cat("Starting Step 4: Normalizing the DGEList and setting the contrasts as defined by Alexa ...\n")  
  
  dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = edgeR_norm)
  
  # Design matrix (no intercept was included since it is a simple two group comparison)
  design_mtx = model.matrix(~ 0 + group)
  
  # Assigns the name of the levels to the design mtx colnames
  colnames(design_mtx) <- levels(group)
  
  # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
  
  if (is.null(contrast_formula)) {
    stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
  }
  
  # the provided colnames must reference the existing design columns
  cn <- colnames(design_mtx)
  
  parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
  
  if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
  
  if (!all(parts %in% cn)) {
    stop(sprintf(
      "Contrast groups not found in design. Provided: %s. Available: %s",
      paste(parts, collapse = ", "),
      paste(cn, collapse = ", ")
    ))
  }
  
  contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
  print(contrasts)
  
  ##### Fit glmQLFit model and test with glmQLFTest the contrasts #####
  cat("Starting Step 5: Fitting voom and lmFit and testing the contrasts ...\n")  
  
  # Fit the quasi-likelihood model
  fit = glmQLFit(dgelist_filtered_norm, design = design_mtx, robust = TRUE)
  
  # Quasi-likelihood F-test
  qlf = glmQLFTest(fit, contrast = contrasts)
  
  # Extract the full results
  de_results = as.data.frame(topTags(qlf, n = Inf))
  
  ##### Create the results table #####
  cat("Starting Step 6: Calculating the metrics and generating the results table ...\n")  
  de_results <- as.data.frame(de_results)
  de_results$X <- rownames(de_results)
  rownames(de_results) <- NULL
  
  # Ensure we only evaluate on tested IDs that appear in results
  # (should match, but we keep this safe)
  tested_ids_try2 <- de_results$X
  truth_tested2 <- truth_set_tested %>% dplyr::filter(X %in% tested_ids_try2)
  truth_tested2 <- truth_tested2 %>% dplyr::distinct(X, .keep_all = TRUE)
  
  # Extract the true DE (TP) and non DE (TN)
  DE_truth <- truth_tested2 %>% dplyr::filter(DE.ind == TRUE)
  DE_notTruth <- truth_tested2 %>% dplyr::filter(DE.ind == FALSE)
  
  # Helper for safe division
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  
  #ROC-AUC, AUPRC
  scored <- de_results %>% 
    dplyr::inner_join(truth_tested2, by = "X") %>%
    dplyr::filter(!is.na(PValue))
  
  scores <- 1 - scored$PValue
  truth_label <- as.integer(scored$DE.ind)
  
  roc_obj     <- pROC::roc(response  = truth_label,
                           predictor = scores,
                           quiet     = TRUE)
  gene_roc_auc <- as.numeric(pROC::auc(roc_obj))
  
  pr_obj      <- PRROC::pr.curve(scores.class0 = scores[truth_label == 1],  # truly DE
                                 scores.class1 = scores[truth_label == 0],  # not DE
                                 curve         = FALSE)
  gene_pr_auc <- pr_obj$auc.integral
  
  # Calculate metrics for every specified threshold
  out <- lapply(thresholds, function(th) {
    
    called <- de_results %>% filter(.data[[threshold_col]] <= th)
    not_called <- de_results %>% filter(.data[[threshold_col]] > th)
    
    TP <- nrow(inner_join(called,     DE_truth,  by = "X"))
    FP <- nrow(inner_join(called,     DE_notTruth, by = "X"))
    TN <- nrow(inner_join(not_called, DE_notTruth, by = "X"))
    FN <- nrow(inner_join(not_called, DE_truth,  by = "X"))
    
    precision <- safe_div(TP, TP + FP)
    recall    <- safe_div(TP, TP + FN)  # = sensitivity = TPR
    specificity <- safe_div(TN, TN + FP)
    fpr <- safe_div(FP, FP + TN)
    fdr <- safe_div(FP, TP + FP)
    f1  <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                  NA_real_,
                  2 * precision * recall / (precision + recall))
    
    data.frame(
      rep_id = rep_id,
      threshold_type = threshold_col,
      threshold = th,
      
      initial_genes_truth = n_genes_initial_truth_set,
      initial_genes_counts = n_genes_initial_counts,
      tested_genes = n_circRNAs_tested,
      pct_retained = pct_retained,
      
      TP = TP, FP = FP, TN = TN, FN = FN,
      
      sensitivity = recall * 100,
      specificity = specificity * 100,
      precision = precision,
      recall = recall,
      gene_roc_auc = gene_roc_auc,
      gene_pr_auc  = gene_pr_auc,
      F1 = f1,
      FDR = fdr * 100,
      FPR = fpr * 100,
      
      n_called = nrow(called)
    )
  }) %>% bind_rows()
  
  runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  out$runtime_sec <- runtime_sec
  cat(sprintf(" Completed in %.2f seconds\n", out$runtime_sec))
  
  list(metrics = out, de_results = de_results)
}

##############
### DESeq2 ###
##############

## ZeroSet ##
# Functions #
run_DESeq2_betaprior_DE0 <- function(bsj_counts_path,
                                     metadata_path,
                                     pvalue_thresholds = c(0.01, 0.05, 0.10),
                                     rep_id = NA,
                                     contrast_formula = NULL) {
  set.seed(42)
  # Start timer
  start_time = Sys.time()
  tryCatch({
    
    # Step 1: Load the data and the metadata
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    if("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # Print the initial gene number
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    
    # Make DESeq2 compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # Step 3: Create the DESeq2 object
    cat("Starting Step 3: Creating DESeqDataSet...\n")
    dds <- DESeqDataSetFromMatrix(countData = bsj_fc,
                                  colData = bsj_fc_md,
                                  design = ~Group)
    
    # Step 4: Relevel to compare cancer vs healthy
    # dds$Group <- relevel(dds$Group, ref = "Breast.Cancer.Tissue")
    
    # Step 5: Run DESeq2
    cat("Starting Step 5: Running DESeq2...\n")
    dds <- DESeq(dds,
                 betaPrior = TRUE)
    
    
    # Extract number of genes after DESeq2 filtering**
    # DESeq2 automatically filters out genes with all zeros or very low counts
    # Genes that pass filtering are those with results
    g1 <- strsplit(contrast_formula, "-")[[1]][1]
    g2 <- strsplit(contrast_formula, "-")[[1]][2]
    res_cancer_healthy <- results(dds, 
                                  contrast = c("Group", 
                                               g1,
                                               g2))
    
    res_cancer_healthy.df <- as.data.frame(res_cancer_healthy)
    
    # Remove genes with NA pvalues (these were filtered by DESeq2)
    n_genes_tested <- sum(!is.na(res_cancer_healthy.df$pvalue))
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 6: Calculate the metrics at each Pvalue threshold
    cat("Starting Step 6: Calculating metrics at each Pvalue threshold...\n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    sig_genes_list <- list()
    # For each p-value threshold
    for (pval in pvalue_thresholds) {
      
      # Calculate the sum of rows that fulfill the pval condition
      n_sig = sum(res_cancer_healthy.df$pvalue < pval, na.rm = TRUE)
      
      # False positive rate
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
      
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      
      # Get genes significant at this threshold
      sig_idx <- res_cancer_healthy.df$pvalue < pval
      sig_genes <- rownames(res_cancer_healthy.df)[sig_idx]
      
      # Store as comma-separated string
      sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval = round(median(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval = round(min(res_cancer_healthy.df$pvalue, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))
    
    return(list(
      metrics = metrics, 
      de_results = res_cancer_healthy.df,
      sig_genes = sig_genes_list
    ))
    
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(list(
      metrics = error_row,
      de_results = NULL,
      sig_genes = list()
    ))
  })
}

run_DESeq2_LRT_DE0 <- function(bsj_counts_path,
                               metadata_path,
                               pvalue_thresholds = c(0.01, 0.05, 0.10),
                               rep_id = NA,
                               contrast_formula = NULL) {
  set.seed(42)
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    if("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # Print the initial gene number
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    
    # Make DESeq2 compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # Step 3: Create the DESeq2 object
    cat("Starting Step 3: Creating DESeqDataSet...\n")
    dds <- DESeqDataSetFromMatrix(countData = bsj_fc,
                                  colData = bsj_fc_md,
                                  design = ~Group)
    
    # Step 4: Relevel to compare cancer vs healthy
    # dds$Group <- relevel(dds$Group, ref = "Breast.Cancer.Tissue")
    
    # Step 5: Run DESeq2
    cat("Starting Step 5: Running DESeq2...\n")
    dds <- DESeq(dds,
                 test = "LRT",
                 reduced = ~ 1)
    
    
    # Extract number of genes after DESeq2 filtering**
    # DESeq2 automatically filters out genes with all zeros or very low counts
    # Genes that pass filtering are those with results
    
    g1 <- strsplit(contrast_formula, "-")[[1]][1]
    g2 <- strsplit(contrast_formula, "-")[[1]][2]
    
    res_cancer_healthy <- results(dds, 
                                  contrast = c("Group", 
                                               g1,
                                               g2))
    
    res_cancer_healthy.df <- as.data.frame(res_cancer_healthy)
    
    # Remove genes with NA pvalues (these were filtered by DESeq2)
    n_genes_tested <- sum(!is.na(res_cancer_healthy.df$pvalue))
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 6: Calculate the metrics at each Pvalue threshold
    cat("Starting Step 6: Calculating metrics at each Pvalue threshold...\n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    sig_genes_list <- list()
    # For each p-value threshold
    for (pval in pvalue_thresholds) {
      
      # Calculate the sum of rows that fulfill the pval condition
      n_sig = sum(res_cancer_healthy.df$pvalue < pval, na.rm = TRUE)
      
      # False positive rate
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      
      # Get genes significant at this threshold
      sig_idx <- res_cancer_healthy.df$pvalue < pval
      sig_genes <- rownames(res_cancer_healthy.df)[sig_idx]
      
      # Store as comma-separated string
      sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval = round(median(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval = round(min(res_cancer_healthy.df$pvalue, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))
    
    return(list(
      metrics = metrics, 
      de_results = res_cancer_healthy.df,
      sig_genes = sig_genes_list
    ))
    
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(list(
      metrics = error_row,
      de_results = NULL,
      sig_genes = list()
    ))
  })
}

run_DESeq2_wald_DE0 <- function(bsj_counts_path,
                                metadata_path,
                                pvalue_thresholds = c(0.01, 0.05, 0.10),
                                rep_id = NA,
                                contrast_formula = NULL) {
  set.seed(42)
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    if("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # Print the initial gene number
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    
    # Make DESeq2 compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # Step 3: Create the DESeq2 object
    cat("Starting Step 3: Creating DESeqDataSet...\n")
    dds <- DESeqDataSetFromMatrix(countData = bsj_fc,
                                  colData = bsj_fc_md,
                                  design = ~Group)
    
    # Step 4: Relevel to compare cancer vs healthy
    # dds$Group <- relevel(dds$Group, ref = "Breast.Cancer.Tissue")
    
    # Step 5: Run DESeq2
    cat("Starting Step 5: Running DESeq2...\n")
    dds <- DESeq(dds)
    
    # Extract number of genes after DESeq2 filtering**
    # DESeq2 automatically filters out genes with all zeros or very low counts
    # Genes that pass filtering are those with results
    g1 <- strsplit(contrast_formula, "-")[[1]][1]
    g2 <- strsplit(contrast_formula, "-")[[1]][2]
    res_cancer_healthy <- results(dds, 
                                  contrast = c("Group", 
                                               g1,
                                               g2))
    
    res_cancer_healthy.df <- as.data.frame(res_cancer_healthy)
    print(res_cancer_healthy.df)
    
    # Remove genes with NA pvalues (these were filtered by DESeq2)
    n_genes_tested <- sum(!is.na(res_cancer_healthy.df$pvalue))
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 6: Calculate the metrics at each Pvalue threshold
    cat("Starting Step 6: Calculating metrics at each Pvalue threshold...\n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    sig_genes_list <- list()
    # For each p-value threshold
    for (pval in pvalue_thresholds) {
      
      # Calculate the sum of rows that fulfill the pval condition
      n_sig = sum(res_cancer_healthy.df$pvalue < pval, na.rm = TRUE)
      
      # False positive rate
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
      
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      
      # Get genes significant at this threshold
      sig_idx <- res_cancer_healthy.df$pvalue < pval
      sig_genes <- rownames(res_cancer_healthy.df)[sig_idx]
      
      # Store as comma-separated string
      sig_genes_list[[paste0("sig_genes_", pval_label)]] <- paste(sig_genes, collapse = ", ")
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval = round(median(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval = round(min(res_cancer_healthy.df$pvalue, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))
    
    return(list(
      metrics = metrics, 
      de_results = res_cancer_healthy.df,
      sig_genes = sig_genes_list
    ))
    
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(list(
      metrics = error_row,
      de_results = NULL,
      sig_genes = list()
    ))
  })
}

# Iterators #
DESeq2_betaprior_DE0 <- function(data_dir,
                                 dataset_name,
                                 file_pattern,
                                 pvalue_thresholds = c(0.01, 0.05, 0.10),
                                 output_dir = NULL,
                                 contrast_formula = NULL,
                                 filter = NULL) {
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    # Run analysis
    result <- run_DESeq2_betaprior_DE0(
      bsj_counts_path = counts_path,
      metadata_path = design_path,
      pvalue_thresholds = pvalue_thresholds,
      rep_id = rep_numbers[i],
      contrast_formula = contrast_formula
    )
    
    # Store result
    results_list[[i]] <- result$metrics
    sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), result$sig_genes)
    # Track success/failure
    if (!is.null(result) && !is.na(result$metrics$n_genes_tested[1])) {
      n_success <- n_success + 1
    } else {
      n_failed <- n_failed + 1
    }
  }
  
  #### Step 5: Combine results ####
  cat("\nCombining results...\n")
  results_df <- do.call(rbind, results_list)
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      .before = 1,
      norm_method = "Median of Ratios"
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_%s_DE0_DESeq2_BetaPrior_summary.csv", dataset_name, filter)
  output_path <- file.path(output_dir, output_filename)
  
  write.csv(results_df, output_path, row.names = FALSE)
  sig_genes_filename <- sprintf("%s_%s_DE0_DESeq2_BetaPrior_sig_genes.csv", 
                                dataset_name, filter)
  write.csv(sig_genes_df, file.path(output_dir, sig_genes_filename), row.names = FALSE)
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(results_df)
}

analyze_DESeq2_LRT_DE0 <- function(data_dir,
                                   dataset_name,
                                   file_pattern,
                                   pvalue_thresholds = c(0.01, 0.05, 0.10),
                                   output_dir = NULL,
                                   contrast_formula = NULL,
                                   filter=NULL) {
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    # Run analysis
    result <- run_DESeq2_LRT_DE0(
      bsj_counts_path = counts_path,
      metadata_path = design_path,
      pvalue_thresholds = pvalue_thresholds,
      rep_id = rep_numbers[i],
      contrast_formula = contrast_formula
    )
    
    # Store result
    results_list[[i]] <- result$metrics
    sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), result$sig_genes)
    # Track success/failure
    if (!is.null(result) && !is.na(result$metrics$n_genes_tested[1])) {
      n_success <- n_success + 1
    } else {
      n_failed <- n_failed + 1
    }
  }
  
  #### Step 5: Combine results ####
  cat("\nCombining results...\n")
  results_df <- do.call(rbind, results_list)
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      .before = 1,
      norm_method = "Median of Ratios"
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_%s_DE0_DESeq2_LRT_summary.csv", dataset_name, filter)
  output_path <- file.path(output_dir, output_filename)
  
  write.csv(results_df, output_path, row.names = FALSE)
  sig_genes_filename <- sprintf("%s_%s_DE0_DESeq2_LRT_sig_genes.csv", 
                                dataset_name, filter)
  write.csv(sig_genes_df, file.path(output_dir, sig_genes_filename), row.names = FALSE)
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(results_df)
}

analyze_deseq2_wald_null <- function(data_dir,
                                     dataset_name,
                                     file_pattern,
                                     pvalue_thresholds = c(0.01, 0.05, 0.10),
                                     output_dir = NULL,
                                     contrast_formula = NULL,
                                     filter=NULL) {
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  sig_genes_list <- vector("list", length = length(all_files))
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    # Run analysis
    result <- run_DESeq2_wald_DE0(
      bsj_counts_path = counts_path,
      metadata_path = design_path,
      pvalue_thresholds = pvalue_thresholds,
      rep_id = rep_numbers[i],
      contrast_formula = contrast_formula
    )
    
    # Store result
    results_list[[i]] <- result$metrics
    sig_genes_list[[i]] <- c(rep_id = as.character(rep_numbers[i]), result$sig_genes)
    # Track success/failure
    if (!is.null(result) && !is.na(result$metrics$n_genes_tested[1])) {
      n_success <- n_success + 1
    } else {
      n_failed <- n_failed + 1
    }
  }
  
  #### Step 5: Combine results ####
  cat("\nCombining results...\n")
  results_df <- do.call(rbind, results_list)
  sig_genes_df <- do.call(rbind, lapply(sig_genes_list, function(x) {
    df <- as.data.frame(t(x), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)
    return(df)
  }))
  sig_genes_df$rep_id <- as.numeric(sig_genes_df$rep_id)
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      .before = 1,
      norm_method = "Median of Ratios"
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_%s_DE0_DESeq2_WaldTest_summary.csv", dataset_name, filter)
  output_path <- file.path(output_dir, output_filename)
  
  write.csv(results_df, output_path, row.names = FALSE)
  sig_genes_filename <- sprintf("%s_%s_DE0_DESeq2_WaldTest_sig_genes.csv", 
                                dataset_name, filter)
  write.csv(sig_genes_df, file.path(output_dir, sig_genes_filename), row.names = FALSE)
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(results_df)
}

## Signal Set ##
# Functions #
run_DESeq2_WaldTest_with_truth <- function(bsj_counts_path,
                                           metadata_path,
                                           truth_set_path,
                                           pvalue_thresholds = c(0.01, 0.05, 0.10),
                                           threshold_col = c("padj"),
                                           rep_id = NA,
                                           contrast_formula = NULL,
                                           quiet = TRUE) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  tryCatch({
    
    # Step 1: Load counts
    cat("Step 1: Loading counts...\n")
    bsj_fc <- read.csv(bsj_counts_path)
    
    if ("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # DESeq2 wants integer counts
    bsj_fc <- as.matrix(bsj_fc)
    storage.mode(bsj_fc) <- "integer"
    
    n_genes_initial <- nrow(bsj_fc)
    
    # Step 2: Load metadata
    cat("Step 2: Loading metadata...\n")
    bsj_fc_md <- read.csv(metadata_path)
    
    # Ensure sample_id exists (many of your design files may not include it)
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      bsj_fc_md$sample_id <- colnames(bsj_fc)
    }
    
    # Reorder metadata to match counts
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Make Group DESeq2-compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # AUTO-DETECT contrast if not provided
    if (is.null(contrast_formula)) {
      group_levels <- levels(bsj_fc_md$Group)
      if (length(group_levels) != 2) {
        stop("Expected exactly 2 group levels, found: ", paste(group_levels, collapse = ", "))
      }
      # Use alphabetical order (DESeq2 default)
      contrast <- c("Group", group_levels[2], group_levels[1])
      cat(sprintf("Auto-detected contrast: %s vs %s\n", group_levels[2], group_levels[1]))
    } else {
      # Parse contrast_formula like "Breast.Cancer.Tissue-Normal.Adjacent.Tissue"
      groups <- trimws(unlist(strsplit(contrast_formula, "-")))
      groups <- make.names(groups)  # Make DESeq2-compatible
      contrast <- c("Group", groups[1], groups[2])
      cat(sprintf("Using provided contrast: %s vs %s\n", groups[1], groups[2]))
    }
    
    # Step 3: Load truth set
    cat("Step 3: Loading truth set...\n")
    truth_set <- read.csv(truth_set_path)
    if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
      stop("Truth file must contain columns: X, DE.ind")
    }
    
    truth_set$DE.ind <- as.logical(truth_set$DE.ind)
    
    # Important: removes duplicates
    truth_set <- truth_set %>% distinct(X, .keep_all = TRUE)
    
    n_truth_initial <- length(unique(truth_set$X))
    
    # Step 4: DESeq2 object
    cat("Step 4: Creating DESeqDataSet...\n")
    
    dds <- DESeqDataSetFromMatrix(
      countData = bsj_fc,
      colData = bsj_fc_md,
      design = ~ Group
    )
    
    
    # Step 5: Run DESeq2
    cat("Step 5: Running DESeq2...\n")
    
    dds <- DESeq(dds)
    
    # Step 6: Results
    cat("Step 6: Extracting results...\n")
    res <- results(dds, contrast = contrast)
    
    res_df <- as.data.frame(res)
    
    # Tested universe: non-NA in the chosen threshold column
    res_df <- as.data.frame(res) %>% rownames_to_column("X")
    tested_df <- res_df %>% filter(!is.na(.data[[threshold_col]]))
    
    # these are the tested circRNAs for DE
    tested_ids <- tested_df$X
    
    # Genes that were tested in the DE analysis by DESeq2
    n_genes_tested <- length(tested_ids)
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf("Retained %d/%d genes (%.1f%%) in tested set\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Restrict truth to tested universe
    truth_tested <- truth_set %>% filter(X %in% tested_ids)
    
    DE_truth     <- truth_tested %>% filter(DE.ind == TRUE)
    nonDE_truth  <- truth_tested %>% filter(DE.ind == FALSE)
    
    # Helpers
    safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
    
    #ROC-AUC, AUPRC
    scored <- tested_df %>%
      dplyr::inner_join(truth_tested, by = "X") %>%
      dplyr::filter(!is.na(pvalue))
    scores      <- 1 - scored$pvalue        
    truth_label <- as.integer(scored$DE.ind)
    
    roc_obj      <- pROC::roc(response  = truth_label,
                              predictor = scores,
                              quiet     = TRUE)
    gene_roc_auc <- as.numeric(pROC::auc(roc_obj))
    
    pr_obj       <- PRROC::pr.curve(scores.class0 = scores[truth_label == 1],
                                    scores.class1 = scores[truth_label == 0],
                                    curve         = FALSE)
    gene_pr_auc  <- pr_obj$auc.integral
    
    
    # Store metrics per threshold (tidy format: one row per threshold)
    metrics <- lapply(pvalue_thresholds, function(th) {
      
      called <- tested_df %>% filter(.data[[threshold_col]] < th)
      not_called <- tested_df %>% filter(!(X %in% called$X))
      
      TP <- nrow(inner_join(called,     DE_truth,    by = "X"))
      FP <- nrow(inner_join(called,     nonDE_truth, by = "X"))
      TN <- nrow(inner_join(not_called, nonDE_truth, by = "X"))
      FN <- nrow(inner_join(not_called, DE_truth,    by = "X"))
      
      precision <- safe_div(TP, TP + FP)
      recall    <- safe_div(TP, TP + FN)
      specificity <- safe_div(TN, TN + FP)
      fpr <- safe_div(FP, FP + TN)
      fdr <- safe_div(FP, TP + FP)
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                   NA_real_,
                   2 * precision * recall / (precision + recall))
      
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        
        n_genes_initial = n_genes_initial,
        n_truth_initial = n_truth_initial,
        n_genes_tested = n_genes_tested,
        pct_filtered = pct_filtered,
        
        TP = TP, FP = FP, TN = TN, FN = FN,
        
        sensitivity = recall * 100,
        specificity = specificity * 100,
        precision = precision,
        recall = recall,
        gene_roc_auc = gene_roc_auc,
        gene_pr_auc  = gene_pr_auc,
        F1 = f1,
        FDR = fdr * 100,
        FPR = fpr * 100,
        
        n_called = nrow(called)
      )
    }) %>% bind_rows()
    
    metrics$mean_pval <- round(mean(res_df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval <- round(median(res_df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval <- round(min(res_df$pvalue, na.rm = TRUE), 6)
    
    metrics$runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
    
    return(list(metrics = metrics, de_results = res_df))
    
  }, error = function(e) {
    
    cat(sprintf("ERROR: %s\n", e$message))
    
    # Return a tidy error table with one row per threshold (helps batch binding)
    error_rows <- lapply(pvalue_thresholds, function(th) {
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        n_genes_initial = NA,
        n_truth_initial = NA,
        n_genes_tested = NA,
        pct_filtered = NA,
        TP = NA, FP = NA, TN = NA, FN = NA,
        sensitivity = NA, specificity = NA,
        precision = NA, recall = NA, F1 = NA,
        FDR = NA, FPR = NA,
        n_called = NA,
        mean_pval = NA, median_pval = NA, min_pval = NA,
        runtime_sec = NA,
        gene_roc_auc = NA,
        gene_pr_auc  = NA,
        error_message = e$message
      )
    }) %>% bind_rows()
    
    return(list(metrics = error_rows))
  })
}

run_DESeq2_LRT_with_truth <- function(bsj_counts_path,
                                      metadata_path,
                                      truth_set_path,
                                      pvalue_thresholds = c(0.01, 0.05, 0.10),
                                      threshold_col = c("padj"),
                                      rep_id = NA,
                                      contrast_formula = NULL,
                                      quiet = TRUE) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  tryCatch({
    
    # Step 1: Load counts
    cat("Step 1: Loading counts...\n")
    bsj_fc <- read.csv(bsj_counts_path)
    
    if ("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # DESeq2 wants integer counts
    bsj_fc <- as.matrix(bsj_fc)
    storage.mode(bsj_fc) <- "integer"
    
    n_genes_initial <- nrow(bsj_fc)
    
    # Step 2: Load metadata
    cat("Step 2: Loading metadata...\n")
    bsj_fc_md <- read.csv(metadata_path)
    
    # Ensure sample_id exists (many of your design files may not include it)
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      bsj_fc_md$sample_id <- colnames(bsj_fc)
    }
    
    # Reorder metadata to match counts
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Make Group DESeq2-compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # AUTO-DETECT contrast if not provided
    if (is.null(contrast_formula)) {
      group_levels <- levels(bsj_fc_md$Group)
      if (length(group_levels) != 2) {
        stop("Expected exactly 2 group levels, found: ", paste(group_levels, collapse = ", "))
      }
      # Use alphabetical order (DESeq2 default)
      contrast <- c("Group", group_levels[2], group_levels[1])
      cat(sprintf("Auto-detected contrast: %s vs %s\n", group_levels[2], group_levels[1]))
    } else {
      # Parse contrast_formula like "Breast.Cancer.Tissue-Normal.Adjacent.Tissue"
      groups <- trimws(unlist(strsplit(contrast_formula, "-")))
      groups <- make.names(groups)  # Make DESeq2-compatible
      contrast <- c("Group", groups[1], groups[2])
      cat(sprintf("Using provided contrast: %s vs %s\n", groups[1], groups[2]))
    }
    
    # Step 3: Load truth set
    cat("Step 3: Loading truth set...\n")
    truth_set <- read.csv(truth_set_path)
    if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
      stop("Truth file must contain columns: X, DE.ind")
    }
    
    truth_set$DE.ind <- as.logical(truth_set$DE.ind)
    
    # Important: removes duplicates
    truth_set <- truth_set %>% distinct(X, .keep_all = TRUE)
    
    n_truth_initial <- length(unique(truth_set$X))
    
    # Step 4: DESeq2 object
    cat("Step 4: Creating DESeqDataSet...\n")
    
    dds <- DESeqDataSetFromMatrix(
      countData = bsj_fc,
      colData = bsj_fc_md,
      design = ~ Group
    )
    
    
    # Step 5: Run DESeq2
    cat("Step 5: Running DESeq2...\n")
    
    dds <- DESeq(dds,
                 test = "LRT",
                 reduced = ~ 1)
    
    
    # Step 6: Results
    cat("Step 6: Extracting results...\n")
    res <- results(dds,
                   contrast = contrast)
    
    res_df <- as.data.frame(res)
    
    # Tested universe: non-NA in the chosen threshold column
    res_df <- as.data.frame(res) %>% rownames_to_column("X")
    tested_df <- res_df %>% filter(!is.na(.data[[threshold_col]]))
    
    # these are the tested circRNAs for DE
    tested_ids <- tested_df$X
    
    # Genes that were tested in the DE analysis by DESeq2
    n_genes_tested <- length(tested_ids)
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf("Retained %d/%d genes (%.1f%%) in tested set\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Restrict truth to tested universe
    truth_tested <- truth_set %>% filter(X %in% tested_ids)
    
    DE_truth     <- truth_tested %>% filter(DE.ind == TRUE)
    nonDE_truth  <- truth_tested %>% filter(DE.ind == FALSE)
    
    # Helpers
    safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
    
    #ROC-AUC, AUPRC
    scored <- tested_df %>%
      dplyr::inner_join(truth_tested, by = "X") %>%
      dplyr::filter(!is.na(pvalue))
    scores      <- 1 - scored$pvalue        
    truth_label <- as.integer(scored$DE.ind)
    
    roc_obj      <- pROC::roc(response  = truth_label,
                              predictor = scores,
                              quiet     = TRUE)
    gene_roc_auc <- as.numeric(pROC::auc(roc_obj))
    
    pr_obj       <- PRROC::pr.curve(scores.class0 = scores[truth_label == 1],
                                    scores.class1 = scores[truth_label == 0],
                                    curve         = FALSE)
    gene_pr_auc  <- pr_obj$auc.integral
    
    
    # Store metrics per threshold (tidy format: one row per threshold)
    metrics <- lapply(pvalue_thresholds, function(th) {
      
      called <- tested_df %>% filter(.data[[threshold_col]] < th)
      not_called <- tested_df %>% filter(!(X %in% called$X))
      
      TP <- nrow(inner_join(called,     DE_truth,    by = "X"))
      FP <- nrow(inner_join(called,     nonDE_truth, by = "X"))
      TN <- nrow(inner_join(not_called, nonDE_truth, by = "X"))
      FN <- nrow(inner_join(not_called, DE_truth,    by = "X"))
      
      precision <- safe_div(TP, TP + FP)
      recall    <- safe_div(TP, TP + FN)
      specificity <- safe_div(TN, TN + FP)
      fpr <- safe_div(FP, FP + TN)
      fdr <- safe_div(FP, TP + FP)
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                   NA_real_,
                   2 * precision * recall / (precision + recall))
      
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        
        n_genes_initial = n_genes_initial,
        n_truth_initial = n_truth_initial,
        n_genes_tested = n_genes_tested,
        pct_filtered = pct_filtered,
        
        TP = TP, FP = FP, TN = TN, FN = FN,
        
        sensitivity = recall * 100,
        specificity = specificity * 100,
        precision = precision,
        recall = recall,
        gene_roc_auc = gene_roc_auc,
        gene_pr_auc  = gene_pr_auc,
        F1 = f1,
        FDR = fdr * 100,
        FPR = fpr * 100,
        
        n_called = nrow(called)
      )
    }) %>% bind_rows()
    
    metrics$mean_pval <- round(mean(res_df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval <- round(median(res_df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval <- round(min(res_df$pvalue, na.rm = TRUE), 6)
    
    metrics$runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
    
    return(list(metrics = metrics, de_results = res_df))
    
  }, error = function(e) {
    
    cat(sprintf("ERROR: %s\n", e$message))
    
    # Return a tidy error table with one row per threshold (helps batch binding)
    error_rows <- lapply(pvalue_thresholds, function(th) {
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        n_genes_initial = NA,
        n_truth_initial = NA,
        n_genes_tested = NA,
        pct_filtered = NA,
        TP = NA, FP = NA, TN = NA, FN = NA,
        sensitivity = NA, specificity = NA,
        precision = NA, recall = NA, F1 = NA,
        FDR = NA, FPR = NA,
        n_called = NA,
        mean_pval = NA, median_pval = NA, min_pval = NA,
        runtime_sec = NA,
        gene_roc_auc = NA,
        gene_pr_auc  = NA,
        error_message = e$message
      )
    }) %>% bind_rows()
    
    return(list(metrics = error_rows))
  })
}

run_DESeq2_BetaPriorT_with_truth <- function(bsj_counts_path,
                                             metadata_path,
                                             truth_set_path,
                                             pvalue_thresholds = c(0.01, 0.05, 0.10),
                                             threshold_col = c("padj"),
                                             rep_id = NA,
                                             contrast_formula = NULL,
                                             quiet = TRUE) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  tryCatch({
    
    # Step 1: Load counts
    cat("Step 1: Loading counts...\n")
    bsj_fc <- read.csv(bsj_counts_path)
    
    if ("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # DESeq2 wants integer counts
    bsj_fc <- as.matrix(bsj_fc)
    storage.mode(bsj_fc) <- "integer"
    
    n_genes_initial <- nrow(bsj_fc)
    
    # Step 2: Load metadata
    cat("Step 2: Loading metadata...\n")
    bsj_fc_md <- read.csv(metadata_path)
    
    # Ensure sample_id exists (many of your design files may not include it)
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      bsj_fc_md$sample_id <- colnames(bsj_fc)
    }
    
    # Reorder metadata to match counts
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Make Group DESeq2-compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # AUTO-DETECT contrast if not provided
    if (is.null(contrast_formula)) {
      group_levels <- levels(bsj_fc_md$Group)
      if (length(group_levels) != 2) {
        stop("Expected exactly 2 group levels, found: ", paste(group_levels, collapse = ", "))
      }
      # Use alphabetical order (DESeq2 default)
      contrast <- c("Group", group_levels[2], group_levels[1])
      cat(sprintf("Auto-detected contrast: %s vs %s\n", group_levels[2], group_levels[1]))
    } else {
      # Parse contrast_formula like "Breast.Cancer.Tissue-Normal.Adjacent.Tissue"
      groups <- trimws(unlist(strsplit(contrast_formula, "-")))
      groups <- make.names(groups)  # Make DESeq2-compatible
      contrast <- c("Group", groups[1], groups[2])
      cat(sprintf("Using provided contrast: %s vs %s\n", groups[1], groups[2]))
    }
    
    # Step 3: Load truth set
    cat("Step 3: Loading truth set...\n")
    truth_set <- read.csv(truth_set_path)
    if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
      stop("Truth file must contain columns: X, DE.ind")
    }
    
    truth_set$DE.ind <- as.logical(truth_set$DE.ind)
    
    # Important: removes duplicates
    truth_set <- truth_set %>% distinct(X, .keep_all = TRUE)
    
    n_truth_initial <- length(unique(truth_set$X))
    
    # Step 4: DESeq2 object
    cat("Step 4: Creating DESeqDataSet...\n")
    
    dds <- DESeqDataSetFromMatrix(
      countData = bsj_fc,
      colData = bsj_fc_md,
      design = ~ Group
    )
    
    
    # Step 5: Run DESeq2
    cat("Step 5: Running DESeq2...\n")
    
    dds <- DESeq(dds,
                 betaPrior = TRUE)
    
    # Step 6: Results
    cat("Step 6: Extracting results...\n")
    res <- results(dds,
                   contrast = contrast)
    
    res_df <- as.data.frame(res)
    
    # Tested universe: non-NA in the chosen threshold column
    res_df <- as.data.frame(res) %>% rownames_to_column("X")
    tested_df <- res_df %>% filter(!is.na(.data[[threshold_col]]))
    
    # these are the tested circRNAs for DE
    tested_ids <- tested_df$X
    
    # Genes that were tested in the DE analysis by DESeq2
    n_genes_tested <- length(tested_ids)
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf("Retained %d/%d genes (%.1f%%) in tested set\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Restrict truth to tested universe
    truth_tested <- truth_set %>% filter(X %in% tested_ids)
    
    DE_truth     <- truth_tested %>% filter(DE.ind == TRUE)
    nonDE_truth  <- truth_tested %>% filter(DE.ind == FALSE)
    
    # Helpers
    safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
    
    #ROC-AUC, AUPRC
    scored <- tested_df %>%
      dplyr::inner_join(truth_tested, by = "X") %>%
      dplyr::filter(!is.na(pvalue))
    scores      <- 1 - scored$pvalue        
    truth_label <- as.integer(scored$DE.ind)
    
    roc_obj      <- pROC::roc(response  = truth_label,
                              predictor = scores,
                              quiet     = TRUE)
    gene_roc_auc <- as.numeric(pROC::auc(roc_obj))
    
    pr_obj       <- PRROC::pr.curve(scores.class0 = scores[truth_label == 1],
                                    scores.class1 = scores[truth_label == 0],
                                    curve         = FALSE)
    gene_pr_auc  <- pr_obj$auc.integral
    
    # Store metrics per threshold (tidy format: one row per threshold)
    metrics <- lapply(pvalue_thresholds, function(th) {
      
      called <- tested_df %>% filter(.data[[threshold_col]] < th)
      not_called <- tested_df %>% filter(!(X %in% called$X))
      
      TP <- nrow(inner_join(called,     DE_truth,    by = "X"))
      FP <- nrow(inner_join(called,     nonDE_truth, by = "X"))
      TN <- nrow(inner_join(not_called, nonDE_truth, by = "X"))
      FN <- nrow(inner_join(not_called, DE_truth,    by = "X"))
      
      precision <- safe_div(TP, TP + FP)
      recall    <- safe_div(TP, TP + FN)
      specificity <- safe_div(TN, TN + FP)
      fpr <- safe_div(FP, FP + TN)
      fdr <- safe_div(FP, TP + FP)
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                   NA_real_,
                   2 * precision * recall / (precision + recall))
      
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        
        n_genes_initial = n_genes_initial,
        n_truth_initial = n_truth_initial,
        n_genes_tested = n_genes_tested,
        pct_filtered = pct_filtered,
        
        TP = TP, FP = FP, TN = TN, FN = FN,
        
        sensitivity = recall * 100,
        specificity = specificity * 100,
        precision = precision,
        recall = recall,
        F1 = f1,
        gene_roc_auc = gene_roc_auc,
        gene_pr_auc  = gene_pr_auc,
        FDR = fdr * 100,
        FPR = fpr * 100,
        
        n_called = nrow(called)
      )
    }) %>% bind_rows()
    
    metrics$mean_pval <- round(mean(res_df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval <- round(median(res_df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval <- round(min(res_df$pvalue, na.rm = TRUE), 6)
    
    metrics$runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
    
    return(list(metrics = metrics, de_results = res_df))
    
  }, error = function(e) {
    
    cat(sprintf("ERROR: %s\n", e$message))
    
    # Return a tidy error table with one row per threshold (helps batch binding)
    error_rows <- lapply(pvalue_thresholds, function(th) {
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        n_genes_initial = NA,
        n_truth_initial = NA,
        n_genes_tested = NA,
        pct_filtered = NA,
        TP = NA, FP = NA, TN = NA, FN = NA,
        sensitivity = NA, specificity = NA,
        precision = NA, recall = NA, F1 = NA,
        FDR = NA, FPR = NA,
        n_called = NA,
        mean_pval = NA, median_pval = NA, min_pval = NA,
        runtime_sec = NA,
        gene_roc_auc = NA,
        gene_pr_auc  = NA,
        error_message = e$message
      )
    }) %>% bind_rows()
    
    return(list(metrics = error_rows))
  })
}

### Main signal dataset iterator w/o Limma-Voom quantile ###
analyze_simulation_batch_with_truth <- function(data_dir,
                                                dataset_name,
                                                file_pattern,
                                                runner_fun,               # <- pass run_limma_voom_with_truth or run_edgeR_with_truth
                                                norm_method = "TMM",
                                                contrast_formula,
                                                thresholds = c(0.01, 0.05, 0.10),
                                                threshold_col = c("adj.P.Val", "P.Value", "PValue", "FDR", "padj"),
                                                output_dir = NULL,
                                                filter_type = NULL, # needed for output
                                                tool = NULL, # Helps determine stuff
                                                method = NULL # only for DESeq2,
                                                ) {
  set.seed(42)
  threshold_col <- match.arg(threshold_col)
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat(sprintf("Normalization method: %s\n", norm_method))
  cat(sprintf("Threshold column: %s\n", threshold_col))
  cat("========================================\n\n")
  
  if (is.null(output_dir)) output_dir <- data_dir
  if (!endsWith(data_dir, "/")) data_dir <- paste0(data_dir, "/")
  
  ##### Step 1: Find all count files #####
  cat("Discovering simulation files...\n")
  count_files <- list.files(
    data_dir,
    pattern = paste0(file_pattern, ".*_counts\\.csv$"),
    full.names = FALSE
  )
  
  if (length(count_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  cat(sprintf("Found %d count files\n\n", length(count_files)))
  
  ##### Step 2: Extract dataset rep numbers and the count files #####
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", count_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Order the number of the replicates
  ord <- order(rep_numbers)
  cat("Ordering replicates:\n")
  cat(sprintf("  Original order: %s\n", paste(rep_numbers, collapse = ", ")))
  
  count_files <- count_files[ord]
  rep_numbers <- rep_numbers[ord]
  cat(sprintf("  Sorted order: %s\n\n", paste(rep_numbers, collapse = ", ")))
  
  # Extract the design and truth files
  design_files <- gsub("_counts\\.csv$", "_design.csv", count_files)
  truth_files  <- gsub("_counts\\.csv$", "_TP.csv", count_files)
  
  # Check existence of these files in the defined data dir path
  design_exist <- file.exists(paste0(data_dir, design_files))
  truth_exist  <- file.exists(paste0(data_dir, truth_files))
  
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  if (!all(truth_exist)) {
    missing <- truth_files[!truth_exist]
    stop(sprintf("Missing TP truth files: %s", paste(missing, collapse = ", ")))
  }
  
  ##### Step 3: Initiate result storing in a list #####
  results_list <- vector("list", length(count_files))
  de_results_list <- vector("list", length(count_files))
  # as long as the number of count files
  n_total <- length(count_files)
  n_success <- 0
  n_failed <- 0
  
  #### Step 4: Process each replicate dataset #####
  for (i in seq_along(count_files)) {
    rep_id <- rep_numbers[i] # gets updated with each iteration
    cat("\n----------------------------------------\n")
    cat(sprintf("Processing replicate %d/%d (rep_%d)\n", i, n_total, rep_id))
    cat("----------------------------------------\n")
    
    counts_path <- paste0(data_dir, count_files[i])
    design_path <- paste0(data_dir, design_files[i])
    truth_path  <- paste0(data_dir, truth_files[i])
    
    # Sanity Check 1: Print which files are being loaded for this iteration
    cat(sprintf("Files for rep_%d:\n", rep_id))
    cat(sprintf("  Counts:  %s\n", count_files[i]))
    cat(sprintf("  Design:  %s\n", design_files[i]))
    cat(sprintf("  Truth:   %s\n", truth_files[i]))
    
    # Sanity Check 2: Confirm files exist and check file sizes
    cat("\nFile verification:\n")
    cat(sprintf("  Counts exists: %s (size: %d bytes)\n", 
                file.exists(counts_path), 
                ifelse(file.exists(counts_path), file.info(counts_path)$size, 0)))
    cat(sprintf("  Design exists: %s (size: %d bytes)\n", 
                file.exists(design_path), 
                ifelse(file.exists(design_path), file.info(design_path)$size, 0)))
    cat(sprintf("  Truth exists:  %s (size: %d bytes)\n", 
                file.exists(truth_path), 
                ifelse(file.exists(truth_path), file.info(truth_path)$size, 0)))
    
    # Sanity Check 3: Extract rep number from each filename to ensure they match
    counts_rep <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", count_files[i])
    design_rep <- gsub(paste0(".*", file_pattern, "(\\d+)_design\\.csv"), "\\1", design_files[i])
    truth_rep  <- gsub(paste0(".*", file_pattern, "(\\d+)_TP\\.csv"), "\\1", truth_files[i])
    
    cat(sprintf("\nReplicate ID consistency check:\n"))
    cat(sprintf("  Expected rep_id: %d\n", rep_id))
    cat(sprintf("  Counts file rep: %s\n", counts_rep))
    cat(sprintf("  Design file rep: %s\n", design_rep))
    cat(sprintf("  Truth file rep:  %s\n", truth_rep))
    
    if (counts_rep != design_rep || counts_rep != truth_rep) {
      warning(sprintf("MISMATCH: Files for iteration %d do not all correspond to the same replicate!", i))
    } else {
      cat(sprintf(" ATTENTION: All files correspond to rep_%d\n", rep_id))
    }
    
    cat("\nRunning analysis...\n")
    
    # write the results file
    res <- tryCatch({
      if (tool == "DESeq2") {
        runner_fun(
          bsj_counts_path  = counts_path,
          metadata_path    = design_path,
          truth_set_path   = truth_path,
          #contrast_formula = contrast_formula,
          #edgeR_norm       = norm_method,   # runner may ignore if not edgeR; keep for consistency
          #thresholds       = thresholds,
          threshold_col    = threshold_col,
          rep_id           = rep_id
        )
        } else{ 
          runner_fun(
          bsj_counts_path  = counts_path,
          metadata_path    = design_path,
          truth_set_path   = truth_path,
          contrast_formula = contrast_formula,
          edgeR_norm       = norm_method,   # runner may ignore if not edgeR; keep for consistency
          thresholds       = thresholds,
          threshold_col    = threshold_col,
          rep_id           = rep_id
      )}
    }, error = function(e) {
      cat(sprintf(" WARNING: ERROR in rep_%d: %s\n", rep_id, e$message))
      return(list(metrics = data.frame(
        rep_id = rep_id,
        error_message = e$message,
        stringsAsFactors = FALSE
      )))
    })
    
    # Store metrics (expected: data.frame with one row per threshold)
    results_list[[i]] <- res$metrics
    de_results_list[[i]] <- res$de_results
    # Update counters
    if ("error_message" %in% colnames(res$metrics)) {
      n_failed <- n_failed + 1
      cat(sprintf("  Status: FAILED\n"))
    } else {
      n_success <- n_success + 1
      cat(sprintf("  Status: SUCCESS\n"))
    }
  }
  
  ##### Step 5: Combine the results per threshold level ####
  cat("\n========================================\n")
  cat("Combining results...\n")
  results_df <- dplyr::bind_rows(results_list)
  names(de_results_list) <- paste0("rep_", rep_numbers)
  ##### Step 6: Add metadata  #####
  results_df <- results_df %>%
    dplyr::mutate(
      dataset = dataset_name,
      norm_method = norm_method,
      .before = 1
    )
  
  de_results_combined <- dplyr::bind_rows(
    lapply(names(de_results_list), function(rep_name) {
      if (!is.null(de_results_list[[rep_name]])) {
        de_results_list[[rep_name]]$rep_id <- rep_name  # Add rep identifier
        return(de_results_list[[rep_name]])
      }
    })
  )
  ##### Step 7: Save the output file #####
  if (tool =="DESeq2") {
   output_filename <- sprintf("%s_%s_%s_%s_%s_summary.csv", dataset_name, filter_type, "DE0.1", tool, method) 

  } else {
  output_filename <- sprintf("%s_%s_%s_%s_%s_summary.csv", dataset_name, filter_type, "DE0.1", tool, norm_method)
  }
  output_path <- file.path(output_dir, output_filename)
  write.csv(results_df, output_path, row.names = FALSE)
  
  
  de_output_filename <- gsub("_summary.csv", "_DE_results_all_reps.csv", output_filename)
  de_output_path <- file.path(output_dir, de_output_filename)
  write.csv(de_results_combined, de_output_path, row.names = FALSE)

  
  ##### Step 8: Print summary #####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  cat(sprintf("\nDE_Results saved to:\n%s\n", de_output_path))
  # Some quick summaries (only if those columns exist)
  if (n_success > 0 && "tested_genes" %in% colnames(results_df)) {
    # Filter out failed replicates if error_message column exists
    if ("error_message" %in% colnames(results_df)) {
      valid <- results_df %>% dplyr::filter(is.na(error_message))
    } else {
      valid <- results_df
    }
    cat("\n--- Summary Statistics (tested universe) ---\n")
    cat(sprintf("Mean genes tested after filtering: %.0f (±%.0f)\n",
                mean(valid$tested_genes, na.rm = TRUE),
                sd(valid$tested_genes, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  return(list(
    metrics    = results_df,
    de_results = de_results_list
  ))
}

#### Real data ####
## Helpers ##
get_ordered_contrasts <- function(meta_path, contrast_map_entry) {
  meta <- read.csv(meta_path)
  
  # Get actual levels after make.names transformation
  actual_levels <- levels(factor(make.names(meta$Group)))
  
  # Transform the expected labels from contrast_map the same way
  case_label    <- make.names(contrast_map_entry[1])
  control_label <- make.names(contrast_map_entry[2])
  
  # Validate both are present
  if (!case_label %in% actual_levels)
    stop(sprintf("Case label '%s' not found in metadata levels: %s",
                 case_label, paste(actual_levels, collapse = ", ")))
  if (!control_label %in% actual_levels)
    stop(sprintf("Control label '%s' not found in metadata levels: %s",
                 control_label, paste(actual_levels, collapse = ", ")))
  
  list(
    case    = case_label,
    control = control_label,
    levels  = actual_levels
  )
}
save_plots <- function(plots, output_dir, prefix) {
  if (is.null(plots)) return(invisible(NULL))
  for (plot_name in names(plots)) {
    fpath <- file.path(output_dir, paste0(prefix, "_", plot_name, ".png"))
    png(fpath, width = ifelse(plot_name == "voom", 1400, 800), height = 600)
    replayPlot(plots[[plot_name]])
    dev.off()
  }
}
get_csv_files <- function(data_dir, dataset_name, filter_type, de_label) {
  
  # list all CSV files in the folder
  files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # counts and TP files
  counts_files <- files[grepl("_counts_.*_counts\\.csv$", files)]
  tp_files     <- files[grepl("_counts_.*_TP\\.csv$", files)]
  
  # extract replicate numbers from filenames
  extract_rep <- function(fname) {
    m <- regmatches(fname, regexpr("_rep[0-9]+_", fname))
    gsub("_rep", "", gsub("_", "", m))
  }
  
  counts_reps <- sapply(counts_files, extract_rep)
  tp_reps     <- sapply(tp_files, extract_rep)
  
  # match counts and TP by replicate
  paired_files <- lapply(counts_reps, function(rep) {
    counts_f <- counts_files[which(counts_reps == rep)]
    tp_f     <- tp_files[which(tp_reps == rep)]
    
    if (length(counts_f) != 1 || length(tp_f) != 1) {
      stop("Mismatch or missing files for replicate ", rep)
    }
    
    list(
      rep = rep,
      counts = read.csv(counts_f, row.names = 1),
      tp     = read.csv(tp_f, row.names = 1)
    )
  })
  
  names(paired_files) <- paste0("rep", counts_reps)
  paired_files
}

standardise_results <- function(df) {
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df)
  
  # FDR / adjusted p-value
  fdr_col <- intersect(c("FDR", "padj", "adj.P.Val"), colnames(df))[1]
  if (!is.na(fdr_col)) {
    df$fdr <- df[[fdr_col]]
  } else {
    pval_col <- intersect(c("PValue", "pvalue", "P.Value"), colnames(df))[1]
    if (!is.na(pval_col)) {
      df$fdr <- p.adjust(df[[pval_col]], method = "BH")
      df$fdr_source <- "adj.pVal"
      message(sprintf("  No FDR column — computed BH from '%s'", pval_col))
    } else {
      message("  No FDR or p-value column found — skipping")
      return(NULL)
    }
  }
  
  # Raw p-value
  pval_col <- intersect(c("PValue", "pvalue", "P.Value"), colnames(df))[1]
  if (!is.na(pval_col)) df$pval_raw <- df[[pval_col]]
  
  # logFC
  lfc_col <- intersect(c("logFC", "log2FoldChange"), colnames(df))[1]
  if (!is.na(lfc_col)) df$logFC_std <- df[[lfc_col]]
  
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  
  df
}
flatten_results <- function(all_results) {
  out <- list()
  
  for (dataset in names(all_results)) {
    dr <- all_results[[dataset]]
    
    # CiriDE
    for (src in c("ciriFSJ", "STAR")) {
      for (norm in names(dr[["ciriDE"]][[src]])) {
        df <- standardise_results(dr[["ciriDE"]][[src]][[norm]]$results)
        if (!is.null(df))
          out[[length(out)+1]] <- list(
            label   = paste0("CiriDE-", src, "-", norm),
            dataset = dataset, tool = "CiriDE",
            source  = src, mode = "glmLRT",
            norm    = norm, min_cpm = NA, df = df
          )
      }
    }
    
    # edgeR FSJ/STAR
    for (src in c("ciriFSJ", "STAR")) {
      for (norm in names(dr[["edgeR"]][[src]])) {
        df <- standardise_results(dr[["edgeR"]][[src]][[norm]]$results)
        if (!is.null(df))
          out[[length(out)+1]] <- list(
            label   = paste0("edgeR-", src, "-", norm),
            dataset = dataset, tool = "edgeR",
            source  = src, mode = "glmQLF",
            norm    = norm, min_cpm = NA, df = df
          )
      }
    }
    # edgeR BSJonly
    for (norm in names(dr[["edgeR"]][["BSJonly"]])) {
      for (mc in names(dr[["edgeR"]][["BSJonly"]][[norm]])) {
        df <- standardise_results(dr[["edgeR"]][["BSJonly"]][[norm]][[mc]]$results)
        if (!is.null(df))
          out[[length(out)+1]] <- list(
            label   = paste0("edgeR-BSJonly-", norm, "-cpm", mc),
            dataset = dataset, tool = "edgeR",
            source  = "BSJonly", mode = "glmQLF",
            norm    = norm, min_cpm = mc, df = df
          )
      }
    }
    
    # limma FSJ/STAR
    for (src in c("ciriFSJ", "STAR")) {
      for (mode in names(dr[["LV"]][[src]])) {
        for (norm in names(dr[["LV"]][[src]][[mode]])) {
          df <- standardise_results(dr[["LV"]][[src]][[mode]][[norm]]$results)
          if (!is.null(df))
            out[[length(out)+1]] <- list(
              label   = paste0("limma-", src, "-", mode, "-", norm),
              dataset = dataset, tool = "limma",
              source  = src, mode = mode,
              norm    = norm, min_cpm = NA, df = df
            )
        }
      }
    }
    # limma BSJonly
    for (mode in names(dr[["LV"]][["BSJonly"]])) {
      for (norm in names(dr[["LV"]][["BSJonly"]][[mode]])) {
        for (mc in names(dr[["LV"]][["BSJonly"]][[mode]][[norm]])) {
          df <- standardise_results(dr[["LV"]][["BSJonly"]][[mode]][[norm]][[mc]]$results)
          if (!is.null(df))
            out[[length(out)+1]] <- list(
              label   = paste0("limma-BSJonly-", mode, "-", norm, "-cpm", mc),
              dataset = dataset, tool = "limma",
              source  = "BSJonly", mode = mode,
              norm    = norm, min_cpm = mc, df = df
            )
        }
      }
    }
    
    # DESeq2
    for (src in c("ciriFSJ", "STAR", "BSJonly")) {
      for (mode in names(dr[["DESeq2"]][[src]])) {
        df <- standardise_results(dr[["DESeq2"]][[src]][[mode]]$results)
        if (!is.null(df))
          out[[length(out)+1]] <- list(
            label   = paste0("DESeq2-", src, "-", mode),
            dataset = dataset, tool = "DESeq2",
            source  = src, mode = mode,
            norm    = "MedOfRatios", min_cpm = NA, df = df
          )
      }
    }
  }
  out
}
build_summary_table <- function(flat, fdr_thresholds = c(0.01, 0.05, 0.10)) {
  rows <- lapply(flat, function(x) {
    lapply(fdr_thresholds, function(th) {
      df  <- x$df
      sig <- df[!is.na(df$fdr) & df$fdr < th, ]
      
      up       <- if ("logFC_std" %in% colnames(sig)) sum(sig$logFC_std > 0, na.rm = TRUE) else NA
      dn       <- if ("logFC_std" %in% colnames(sig)) sum(sig$logFC_std < 0, na.rm = TRUE) else NA
      avg_pval <- if ("pval_raw" %in% colnames(df)) round(mean(df$pval_raw,   na.rm = TRUE), 6) else NA
      med_pval <- if ("pval_raw" %in% colnames(df)) round(median(df$pval_raw, na.rm = TRUE), 6) else NA
      avg_fdr  <- round(mean(df$fdr,   na.rm = TRUE), 6)
      med_fdr  <- round(median(df$fdr, na.rm = TRUE), 6)
      
      data.frame(
        dataset   = x$dataset,
        tool      = x$tool,
        source    = x$source,   # ciriFSJ / STAR / BSJonly
        mode      = x$mode,     # glmLRT / glmQLF / Default / DFRBT / LmFit / WaldTest / LRT / BetaPrior
        norm      = x$norm,
        min_cpm   = x$min_cpm,
        fdr_threshold = th,
        n_sig     = nrow(sig),
        n_up      = up,
        n_down    = dn,
        avg_pval  = avg_pval,
        med_pval  = med_pval,
        avg_fdr   = avg_fdr,
        med_fdr   = med_fdr,
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}
####### BSJE+FSJ Functions #####
run_ciriDE <- function(bsj_path, fsj_path, meta_path, norm_method, contrasts) {
  
  # Prep the data
  bsj_data <- read.csv(bsj_path, sep = '\t', row.names = 1)
  if ("gene_id" %in% colnames(bsj_data)) {
    gene_map <- setNames(bsj_data$gene_id, rownames(bsj_data))
    bsj_data$gene_id <- NULL
  } else {
    gene_map <- NULL
  }
  
  fsj_data <- read.csv(fsj_path, sep = '\t', row.names = 1)
  metadata <- read.csv(meta_path, sep = ',')

  #Prep FSJ data
  fsj_dge <- DGEList(counts = fsj_data, group = metadata$Group)
  fsj_keep <- filterByExpr(fsj_dge)
  fsj_dge <- fsj_dge[fsj_keep, keep.lib.sizes=F]
  fsj_dge <- calcNormFactors(fsj_dge, method=norm_method)
  
  #Create circDGE
  circ_dge <- DGEList(counts = bsj_data,
                      group=metadata$Group,
                      lib.size=fsj_dge$samples[,"lib.size"],
                      norm.factors= fsj_dge$samples[, "norm.factors"])
  
  # Design matrix — faithful to original: intercept + treat
  ct <- get_ordered_contrasts(meta_path, contrasts)
  treat  <- factor(make.names(metadata$Group), levels = c(ct$control, ct$case))
  design <- model.matrix(~ treat)
  
  # Run the DGE
  circ_dge <- estimateDisp(circ_dge, design, robust=T)
  
  dev.new(visible = FALSE)
  plotBCV(circ_dge, main = paste("CiriDE BCV -", norm_method))
  bcv_plot <- recordPlot()
  dev.off()
  
  circ_fit <- glmFit(circ_dge, design)
  circ_lrt <- glmLRT(circ_fit)
  
  results <- circ_lrt$table
  circ_order <- order(circ_lrt$table$PValue)
  results$DE <- as.vector(decideTests(circ_lrt))
  results <- results[circ_order, ]
  results$adj.PVal <- p.adjust(results$PValue, method="fdr")
  
  return(list(results = results, gene_map = gene_map, plots = list(bcv = bcv_plot)))
}

run_edger <- function(bsj_path, fsj_path, meta_path, norm_method, contrasts) {

  # Prep the data
  bsj_data <- read.csv(bsj_path, sep = '\t', row.names = 1)
  if ("gene_id" %in% colnames(bsj_data)) {
    gene_map <- setNames(bsj_data$gene_id, rownames(bsj_data))
    bsj_data$gene_id <- NULL
  } else {
    gene_map <- NULL
  }
  
  fsj_data <- read.csv(fsj_path, sep = '\t', row.names = 1)
  metadata <- read.csv(meta_path, sep = ',')
  
  #Prep FSJ data
  fsj_dge <- DGEList(counts = fsj_data, group = metadata$Group)
  fsj_keep <- filterByExpr(fsj_dge)
  fsj_dge <- fsj_dge[fsj_keep, keep.lib.sizes=F]
  fsj_dge <- calcNormFactors(fsj_dge, method=norm_method)
  
  #Design mtx
  Group <- factor(make.names(metadata$Group))
  design <- model.matrix(~ 0 + Group)
  colnames(design) <- levels(Group)
  
  ct <- get_ordered_contrasts(meta_path, contrasts)
  
  contrast_vec <- makeContrasts(
    contrasts = paste0(ct$case, "-", ct$control),
    levels = design
  )
  
  #Create circDGE
  circ_dge <- DGEList(counts = bsj_data, group=Group, lib.size=fsj_dge$samples[,"lib.size"], norm.factors= fsj_dge$samples[, "norm.factors"])
  circ_dge <- estimateDisp(circ_dge, design, robust = TRUE)

  # Run the DGE
  fit <- glmQLFit(circ_dge, design = design, robust=T)
  
  dev.new(visible = FALSE)
  plotBCV(circ_dge, main = paste("edgeR BCV -", norm_method))
  bcv_plot <- recordPlot()
  dev.off()
  
  dev.new(visible = FALSE)
  plotQLDisp(fit, main = paste("edgeR QL Dispersion -", norm_method))
  ql_plot <- recordPlot()
  dev.off()
  
  qlf <- glmQLFTest(fit, contrast = contrast_vec)
  
  de_results <- as.data.frame(topTags(qlf, n=Inf))
  
  return(list(results = de_results, gene_map = gene_map, plots = list(bcv = bcv_plot, ql_disp = ql_plot)))
}

run_limma <- function(bsj_path, fsj_path, meta_path, norm_method, contrasts, run_mode) {
  
  # --- Load ---
  bsj <- read.csv(bsj_path, sep = "\t", row.names = 1)
  gene_map <- NULL
  if ("gene_id" %in% colnames(bsj)) {
    gene_map <- setNames(bsj$gene_id, rownames(bsj))
    bsj$gene_id <- NULL
  }
  
  fsj  <- read.csv(fsj_path, sep = "\t", row.names = 1)
  meta <- read.csv(meta_path)
  meta$Group <- factor(make.names(meta$Group))
  
  # --- FSJ normalization ---
  fsj_dge <- DGEList(fsj, group = meta$Group)
  fsj_dge <- fsj_dge[filterByExpr(fsj_dge), , keep.lib.sizes = FALSE]
  fsj_dge <- calcNormFactors(fsj_dge, method = norm_method)
  
  # --- BSJ with borrowed scaling ---
  circ_dge <- DGEList(
    counts       = bsj,
    lib.size     = fsj_dge$samples$lib.size,
    norm.factors = fsj_dge$samples$norm.factors
  )
  
  Group  <- factor(make.names(meta$Group))
  design <- model.matrix(~ 0 + Group)
  colnames(design) <- levels(Group)
  
  ct           <- get_ordered_contrasts(meta_path, contrasts)
  contrast_vec <- makeContrasts(
    contrasts = paste0(ct$case, "-", ct$control),
    levels    = design
  )
  
  # --- Voom + fit ---
  fit_result <- switch(
    run_mode,
    
    Default = {
      dev.new(visible = FALSE)
      v <- voom(circ_dge, design, normalize.method = "none", plot = TRUE)
      voom_plot <- recordPlot()
      dev.off()
      
      fit <- lmFit(v, design)
      fit <- contrasts.fit(fit, contrast_vec)
      list(fit = eBayes(fit), plots = list(voom = voom_plot))
    },
    
    DFRBT = {
      dev.new(visible = FALSE)
      v <- voom(circ_dge, design, normalize.method = "none", plot = TRUE)
      voom_plot <- recordPlot()
      dev.off()
      
      fit <- lmFit(v, design, robust = TRUE)
      fit <- contrasts.fit(fit, contrast_vec)
      list(fit = eBayes(fit), plots = list(voom = voom_plot))
    },
    
    LmFit = {
      dev.new(visible = FALSE)
      par(mfrow = c(1, 2))
      fit <- voomLmFit(circ_dge, design = design, sample.weights = TRUE, plot = TRUE)
      voom_plot <- recordPlot()
      dev.off()
      par(mfrow = c(1, 1))
      
      fit <- contrasts.fit(fit, contrast_vec)
      list(fit = eBayes(fit), plots = list(voom = voom_plot))
    },
    
    stop("Invalid limma mode")
  )
  
  res <- topTable(fit_result$fit, number = Inf, sort.by = "P")
  
  return(list(
    results  = res,
    gene_map = gene_map,
    plots    = fit_result$plots,
    meta     = list(tool = "limma-voom", norm_method = norm_method, mode = run_mode)
  ))
}

run_deseq <- function(bsj_path, fsj_path, meta_path, norm_method, contrasts, run_mode) {
  
  # --- Load ---
  bsj <- read.csv(bsj_path, sep = '\t', row.names = 1)
  if ("gene_id" %in% colnames(bsj)) {
    gene_map <- setNames(bsj$gene_id, rownames(bsj))
    bsj$gene_id <- NULL
  } else {
    gene_map <- NULL
  }
  fsj <- read.csv(fsj_path, sep = "\t", row.names = 1)
  meta <- read.csv(meta_path)
  
  meta$Group <- factor(make.names(meta$Group))
  ct <- get_ordered_contrasts(meta_path, contrasts)
  
  # --- FSJ normalization edgeR method ---
  # fsj_dge <- DGEList(fsj, group = meta$Group)
  # fsj_dge <- fsj_dge[filterByExpr(fsj_dge), , keep.lib.sizes = FALSE]
  # fsj_dge <- calcNormFactors(fsj_dge, method = norm_method)
  # 
  # eff_lib <- fsj_dge$samples$lib.size * fsj_dge$samples$norm.factors
  # size_factors <- eff_lib / exp(mean(log(eff_lib)))
  
  
  # --- FSJ deseq object ---
  fsj_dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(fsj)),
    colData = meta,
    design = ~ Group
  )
  
  fsj_dds     <- estimateSizeFactors(fsj_dds)
  size_factors <- sizeFactors(fsj_dds)
  
  # --- BSJ DESeq2 object ---
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(bsj)),
    colData = meta,
    design = ~ Group
  )
  
  sizeFactors(dds) <- size_factors
  
  # --- Mode switch ---
  dds <- switch(
    run_mode,
    WaldTest = DESeq(dds, test = "Wald"),
    LRT      = DESeq(dds, test = "LRT", reduced = ~ 1),
    BetaPrior= DESeq(dds, betaPrior = TRUE),
    stop("Unknown DESeq2 mode")
  )
  
  dev.new(visible = FALSE)
  plotDispEsts(dds, main = paste("DESeq2 FSJ Dispersion -", run_mode))
  disp_plot <- recordPlot()
  dev.off()
  
  res <- results(dds, contrast = c("Group", ct$case, ct$control))
  res_df <- as.data.frame(res)
  
  return(list(results = res_df, gene_map = gene_map, plots = list(dispersion = disp_plot),
              meta = list(tool = "DESeq2", norm_method = "FSJ-MedianOfRatios", mode = run_mode)))
}

####### BSJ only functions for comparison #####
ciriDE_BSJ <- function(bsj_path, meta_path, norm_method, min_cpm, contrasts) {
  # Prep the data
  bsj_data <- read.csv(bsj_path, sep = '\t', row.names = 1)
  if ("gene_id" %in% colnames(bsj_data)) {
    gene_map <- setNames(bsj_data$gene_id, rownames(bsj_data))
    bsj_data$gene_id <- NULL
  } else {
    gene_map <- NULL
  }
  
  metadata <- read.csv(meta_path, sep = ',')
  
  #Design mtx
  Group <- factor(make.names(metadata$Group))
  design <- model.matrix(~ 0 + Group)
  colnames(design) <- levels(Group)
  
  ct <- get_ordered_contrasts(meta_path, contrasts)
  
  contrast_vec <- makeContrasts(
    contrasts = paste0(ct$case, "-", ct$control),
    levels = design
  )
  
  circ_dge <- DGEList(counts = bsj_data, group = Group)
  if (min_cpm == "def") {
    keep <- filterByExpr(circ_dge, group = Group)
  } else {
    keep <- filterByExpr(circ_dge, group = Group, min.count = as.numeric(min_cpm))
  }
  circ_dge <- circ_dge[keep, , keep.lib.sizes = FALSE]
  circ_dge <- calcNormFactors(circ_dge, method = norm_method)
  
  # Run the DGE
  circ_dge <- estimateDisp(circ_dge, design, robust=T)
  
  dev.new(visible = FALSE)
  plotBCV(circ_dge, main = paste("CiriDE BSJ-only BCV -", norm_method))
  bcv_plot <- recordPlot()
  dev.off()
  
  circ_fit <- glmFit(circ_dge, design)
  circ_lrt <- glmLRT(circ_fit, coef = 2)
  
  results <- circ_lrt$table
  circ_order <- order(circ_lrt$table$PValue)
  results$DE <- as.vector(decideTests(circ_lrt))
  results <- results[circ_order, ]
  results$adj.PVal <- p.adjust(results$PValue, method="fdr")
  
  return(list(results = results, gene_map = gene_map, plots = list(bcv = bcv_plot)))
}

edger_bsj <- function(bsj_path, meta_path, norm_method, min_cpm, contrasts) {
  # Prep the data
  bsj_data <- read.csv(bsj_path, sep = '\t', row.names = 1)
  if ("gene_id" %in% colnames(bsj_data)) {
    gene_map <- setNames(bsj_data$gene_id, rownames(bsj_data))
    bsj_data$gene_id <- NULL
  } else {
    gene_map <- NULL
  }

  metadata <- read.csv(meta_path, sep = ',')
  
  #Design mtx
  Group <- factor(make.names(metadata$Group))
  design <- model.matrix(~ 0 + Group)
  colnames(design) <- levels(Group)
  
  ct <- get_ordered_contrasts(meta_path, contrasts)
  
  contrast_vec <- makeContrasts(
    contrasts = paste0(ct$case, "-", ct$control),
    levels = design
  )
  
  dgelist <- DGEList(counts = bsj_data, group = Group)
  if (min_cpm == "def") {
    keep <- filterByExpr(dgelist, group = Group)
  } else {
    keep <- filterByExpr(dgelist, group = Group, min.count = as.numeric(min_cpm))
  }
  
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = norm_method)
  
  dgelist_filtered_norm <- estimateDisp(dgelist_filtered_norm, design, robust = TRUE)
  # Fit the quasi-likelihood model
  fit = glmQLFit(dgelist_filtered_norm, design = design, robust = TRUE)
  
  dev.new(visible = FALSE)
  plotBCV(dgelist_filtered_norm, main = paste("edgeR BSJonly BCV -", norm_method, "cpm", min_cpm))
  bcv_plot <- recordPlot()
  dev.off()
  
  dev.new(visible = FALSE)
  plotQLDisp(fit, main = paste("edgeR BSJonly QL Dispersion -", norm_method, "cpm", min_cpm))
  ql_plot <- recordPlot()
  dev.off()
  
  # Quasi-likelihood F-test
  qlf = glmQLFTest(fit, contrast = contrast_vec)
  # Extract the full results
  de_results = as.data.frame(topTags(qlf, n = Inf))
  
  return(list(results = de_results, gene_map = gene_map, plots = list(bcv = bcv_plot, ql_disp = ql_plot),
              meta = list(tool = "edgeR", norm_method = norm_method)))
}

limma_bsj <- function(bsj_path, meta_path, norm_method, min_cpm, contrasts, run_mode) {
  
  # --- Load ---
  bsj_data <- read.csv(bsj_path, sep = '\t', row.names = 1)
  if ("gene_id" %in% colnames(bsj_data)) {
    gene_map <- setNames(bsj_data$gene_id, rownames(bsj_data))
    bsj_data$gene_id <- NULL
  } else {
    gene_map <- NULL
  }
  
  metadata <- read.csv(meta_path, sep = ',')
  
  # --- Design ---
  Group  <- factor(make.names(metadata$Group))
  design <- model.matrix(~ 0 + Group)
  colnames(design) <- levels(Group)
  
  ct           <- get_ordered_contrasts(meta_path, contrasts)
  contrast_vec <- makeContrasts(
    contrasts = paste0(ct$case, "-", ct$control),
    levels    = design
  )
  
  # --- Filter + normalize BSJ ---
  dgelist <- DGEList(counts = bsj_data, group = Group)
  keep <- if (min_cpm == "def") {
    filterByExpr(dgelist, group = Group)
  } else {
    filterByExpr(dgelist, group = Group, min.count = as.numeric(min_cpm))
  }
  
  dgelist_filtered      <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = norm_method)
  
  # --- Voom + fit ---
  fit_result <- switch(
    run_mode,
    
    Default = {
      dev.new(visible = FALSE)
      v <- voom(dgelist_filtered_norm, design, normalize.method = "none", plot = TRUE)
      voom_plot <- recordPlot()
      dev.off()
      
      fit <- lmFit(v, design)
      fit <- contrasts.fit(fit, contrast_vec)
      list(fit = eBayes(fit), plots = list(voom = voom_plot))
    },
    
    DFRBT = {
      dev.new(visible = FALSE)
      v <- voom(dgelist_filtered_norm, design, normalize.method = "none", plot = TRUE)
      voom_plot <- recordPlot()
      dev.off()
      
      fit <- lmFit(v, design, robust = TRUE)
      fit <- contrasts.fit(fit, contrast_vec)
      list(fit = eBayes(fit), plots = list(voom = voom_plot))
    },
    
    LmFit = {
      dev.new(visible = FALSE)
      par(mfrow = c(1, 2))
      fit <- voomLmFit(dgelist_filtered_norm, design = design, sample.weights = TRUE, plot = TRUE)
      voom_plot <- recordPlot()
      dev.off()
      par(mfrow = c(1, 1))
      
      fit <- contrasts.fit(fit, contrast_vec)
      list(fit = eBayes(fit), plots = list(voom = voom_plot))
    },
    
    stop("Invalid limma mode")
  )
  
  res <- topTable(fit_result$fit, number = Inf, sort.by = "P")
  
  return(list(
    results  = res,
    gene_map = gene_map,
    plots    = fit_result$plots,
    meta     = list(tool = "limma-voom", norm_method = norm_method, mode = run_mode)
  ))
}

deseq_bsj <- function(bsj_path, meta_path, contrasts, run_mode) {
  
  # --- Load ---
  bsj <- read.csv(bsj_path, sep = '\t', row.names = 1)
  if ("gene_id" %in% colnames(bsj)) {
    gene_map <- setNames(bsj$gene_id, rownames(bsj))
    bsj$gene_id <- NULL
  } else {
    gene_map <- NULL
  }
  meta <- read.csv(meta_path)
  
  meta$Group <- factor(make.names(meta$Group))
  ct <- get_ordered_contrasts(meta_path, contrasts)
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(bsj)),
    colData = meta,
    design = ~ Group
  )
  # --- Mode switch ---
  dds <- switch(
    run_mode,
    WaldTest = DESeq(dds, test = "Wald"),
    LRT      = DESeq(dds, test = "LRT", reduced = ~ 1),
    BetaPrior= DESeq(dds, betaPrior = TRUE),
    stop("Unknown DESeq2 mode")
  )
  
  dev.new(visible = FALSE)
  plotDispEsts(dds, main = paste("DESeq2 BSJonly Dispersion -", run_mode))
  disp_plot <- recordPlot()
  dev.off()
  
  res <- results(dds, contrast = c("Group", ct$case, ct$control))
  res_df <- as.data.frame(res)
  
  return(list(results = res_df, gene_map = gene_map, plots = list(dispersion = disp_plot),
              meta = list(tool = "DESeq2", norm_method = "Median of Ratios", mode = run_mode)))
}


###### Datasets: Paths, maps & helper data ######

#### Sim data ####
### Null data ###
## BSJ only ##
# Autofilter #
bsj_NULL_autofilter_data <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/BC-tissue/autofilter_counts/DE_0",
  EBC1 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC1/autofilter_counts/DE_0",
  EBC2 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC2/autofilter_counts/DE_0",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-PBMC/autofilter_counts/DE_0",
  `HCC-tissue` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-tissue/autofilter_counts/DE_0"
)

# MIN 5 filter #
bsj_NULL_min5_data <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/BC-tissue/min5_counts/DE_0",
  EBC1 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC1/min5_counts/DE_0",
  EBC2 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC2/min5_counts/DE_0",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-PBMC/min5_counts/DE_0",
  `HCC-tissue` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-tissue/min5_counts/DE_0"
)

# MIN 1 filter #
bsj_NULL_min1_data <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/BC-tissue/min1_counts/DE_0",
  EBC1 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC1/min1_counts/DE_0",
  EBC2 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC2/min1_counts/DE_0",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-PBMC/min1_counts/DE_0",
  `HCC-tissue` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-tissue/min1_counts/DE_0"
)

### 10p data ###
## BSJ only ##
# Autofilter #
bsj_10P_autofilter_data <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/BC-tissue/autofilter_counts/DE_0.1",
  EBC1 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC1/autofilter_counts/DE_0.1",
  EBC2 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC2/autofilter_counts/DE_0.1",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-PBMC/autofilter_counts/DE_0.1",
  `HCC-tissue` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-tissue/autofilter_counts/DE_0.1"
)

# MIN 5 filter #
bsj_10P_min5_data <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/BC-tissue/min5_counts/DE_0.1",
  EBC1 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC1/min5_counts/DE_0.1",
  EBC2 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC2/min5_counts/DE_0.1",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-PBMC/min5_counts/DE_0.1",
  `HCC-tissue` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-tissue/min5_counts/DE_0.1"
)

# MIN 1 filter #
bsj_10P_min1_data <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/BC-tissue/min1_counts/DE_0.1",
  EBC1 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC1/min1_counts/DE_0.1",
  EBC2 = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/EBC2/min1_counts/DE_0.1",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-PBMC/min1_counts/DE_0.1",
  `HCC-tissue` = "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only/HCC-tissue/min1_counts/DE_0.1"
)


#### REAL data ####
#BSJ paths
bsj_paths <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/PRJNA553624_breast_cancer/circRNA_outs/CIRI3/BC_CIRI-Candidate.BSJ_Matrix",
  EBC1= "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-11/circRNA_outs/CIRI3/EBC1_CIRI-Candidate.BSJ_Matrix",
  EBC2= "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-25/circRNA_outs/CIRI3/EBC2_CIRI-Candidate.BSJ_Matrix",
  `HCC-PBMC`= "/media/meteor/ChunkyBoi/BC-alarm/CircRNA/PBMC/38samples/circRNA_outs/CIRI3/HCC-PBMC_CIRI-Candidate.BSJ_Matrix",
  `HCC-tissue`= "/media/meteor/FatDawg/Benchmark_Paper/PRJNA716508_hcc_study2/circRNA_outs/CIRI3/HCC-tissue_CIRI-Candidate.BSJ_Matrix"
)

#FSJ paths
fsj_paths <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/PRJNA553624_breast_cancer/circRNA_outs/CIRI3/CIRI-Candidate.FSJ_Matrix",
  EBC1= "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-11/circRNA_outs/CIRI3/CIRI-Candidate.FSJ_Matrix",
  EBC2= "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-25/circRNA_outs/CIRI3/CIRI-Candidate.FSJ_Matrix",
  `HCC-PBMC`= "/media/meteor/ChunkyBoi/BC-alarm/CircRNA/PBMC/38samples/circRNA_outs/CIRI3/CIRI-Candidate.FSJ_Matrix",
  `HCC-tissue`= "/media/meteor/FatDawg/Benchmark_Paper/PRJNA716508_hcc_study2/circRNA_outs/CIRI3/CIRI-Candidate.FSJ_Matrix"
)

#STAR paths
star_fc_paths <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/PRJNA553624_breast_cancer/featureCounts_matrix.txt",
  EBC1= "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-11/2025-11-11-LINEAR-featureCounts_matrix.txt",
  EBC2= "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-25/2025-11-25-LINEAR-featureCounts_matrix.txt",
  `HCC-PBMC`= "/media/meteor/ChunkyBoi/BC-alarm/CircRNA/PBMC/38samples/featureCounts_matrix.txt",
  `HCC-tissue`= "/media/meteor/FatDawg/Benchmark_Paper/PRJNA716508_hcc_study2/featureCounts_matrix.txt"
)

#metadata
metadata_paths <- list(
  BC = "/media/meteor/FatDawg/Benchmark_Paper/metadata_online_datasets/PRJNA553624_metadata.csv",
  EBC1 = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-11/11-11-Metadata.csv",
  EBC2 = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2025-11-25/11-25-Metadata.csv",
  `HCC-PBMC`= "/media/meteor/ChunkyBoi/BC-alarm/CircRNA/PBMC/38samples/PBMC-meta.csv",
  `HCC-tissue`= "/media/meteor/FatDawg/Benchmark_Paper/metadata_online_datasets/PRJNA716508_metadata_cleaned.csv"
)

############# Output directories #############
output_dir <- "/media/meteor/FatDawg/Benchmark_Paper/DE_Results"
ciri_fsj_output_dir <- "/media/meteor/FatDawg/Benchmark_Paper/DE_Results/BSJ-FSJ/ciriDE"
star_fsj_output_dir <- "/media/meteor/FatDawg/Benchmark_Paper/DE_Results/BSJ-FSJ/STAR"

############# Normalization methods ###########
norm_methods <- c('TMM', 'TMMwsp')
de_label_map <- c(
  Null = "DE0",
  `10P` = "DE0.1"
)

contrast_map <- c(
  BC = "Breast.Cancer.Tissue-Normal.Adjacent.Tissue",
  EBC1 = "EBC-Healthy",
  EBC2 = "EBC-Healthy",
  `HCC-PBMC` = "hcc-healthy",
  `HCC-tissue` = "HCC.Adjacent-HCC.Tumor"
)

threshold_col_map <- c(
  limma = "adj.P.val",
  edgeR = "adj.P.val",
  DESeq2 = "FDR"
)

tool_map <- c(
  limma = "Limma-voom",
  DESeq2 = "DESeq2",
  edgeR = "edgeR"
)

variants <- c("default", "LmFit", "DFRBT")

contrast_map <- list(
  BC = c("cancer.tissue", "adjacent.normal.tissue"),
  EBC1 = c("Cancer", "Healthy"),
  EBC2 = c("Cancer", "Healthy"),
  `HCC-PBMC` = c("Hepatocellular.carcinoma", "healthy"),
  `HCC-tissue` = c("tumor", "Adjacent.tissue")
)

limma_modes <- c("Default", "DFRBT", "LmFit") # Quantile is not compatible with this method 
deseq_modes <- c("WaldTest", "LRT", "BetaPrior")

###### Run sim data ######
##### LIMMA loop #######
run_limma_voom <- function(
    data_dir,
    dataset_name,
    filter_type,
    de_label,
    norm,
    contrast_formula,
    output_dir,
    variant
) {
  
  file_pattern <- paste0(dataset_name, "_", filter_type, "_counts_", de_label, "_rep")
  is_null <- de_label == "DE0"
  
  default_res <- if (is_null) {
    limma_voom_TMM_default(
      data_dir = data_dir,
      dataset_name = dataset_name,
      file_pattern = file_pattern,
      norm_method = norm,
      contrast_formula = contrast_formula,
      output_dir = output_dir,
      filter = filter_type,
      de_label = de_label,
      variant = "default"
    )
  } else {
    analyze_simulation_batch_with_truth(
      data_dir = data_dir,
      dataset_name = dataset_name,
      file_pattern = file_pattern,
      runner_fun = run_limma_vdefault_with_truth,
      norm_method = norm,
      contrast_formula = contrast_formula,
      threshold_col = "adj.P.Val",
      output_dir = output_dir,
      filter = filter_type,
      tool = paste0("Limma-voom_", variant)
    )
  }
  
  lmfit_res <- if (is_null) {
    limma_voom_TMM_LmFit_weights(
      data_dir = data_dir,
      dataset_name = dataset_name,
      file_pattern = file_pattern,
      norm_method = norm,
      contrast_formula = contrast_formula,
      output_dir = output_dir,
      filter = filter_type,
      de_label = de_label,
      variant = "LmFit"
    )
  } else {
    analyze_simulation_batch_with_truth(
      data_dir = data_dir,
      dataset_name = dataset_name,
      file_pattern = file_pattern,
      runner_fun = run_limma_LmFit_with_truth,
      norm_method = norm,
      contrast_formula = contrast_formula,
      threshold_col = "adj.P.Val",
      output_dir = output_dir,
      filter = filter_type,
      tool = paste0("Limma-voom_", variant)
    )
  }
  
  dfrbt_res <- if (is_null) {
    limma_voom_DFRBT(
      data_dir = data_dir,
      dataset_name = dataset_name,
      file_pattern = file_pattern,
      norm_method = norm,
      contrast_formula = contrast_formula,
      output_dir = output_dir,
      filter = filter_type,
      de_label = de_label,
      variant = "DFRBT"
    )
  } else {
    analyze_simulation_batch_with_truth(
      data_dir = data_dir,
      dataset_name = dataset_name,
      file_pattern = file_pattern,
      runner_fun = run_limma_vDFBRT_with_truth,
      norm_method = norm,
      contrast_formula = contrast_formula,
      threshold_col = "adj.P.Val",
      output_dir = output_dir,
      filter = filter_type,
      tool = paste0("Limma-voom_", variant)
    )
  }
  
  list(
    voom_default = default_res,
    voom_lmfit   = lmfit_res,
    voom_DFRBT   = dfrbt_res
  )
}

run_limma_quantile <- function(
    data_dir,
    dataset_name,
    filter_type,
    de_label,
    contrast_formula,
    output_dir) {
  
  file_pattern <- paste0(dataset_name, "_", filter_type, "_counts_", de_label, "_rep")
  
  if (de_label == "DE0") {
    limma_voom_quantile(
      data_dir = data_dir,
      dataset_name = dataset_name,
      file_pattern = file_pattern,
      contrast_formula = contrast_formula,
      output_dir = output_dir,
      filter = filter_type,
      de_label = de_label
    )
  } else {
    lv_quantile_signal(
      data_dir = data_dir,
      dataset_name = dataset_name,
      runner_fun = run_limma_vquantile_with_truth,
      file_pattern = file_pattern,
      contrast_formula = contrast_formula,
      output_dir = output_dir,
      filter_type = filter_type
    )
  }
}

limma_start <- Sys.time()
limma_results <- list()

for (de_type in names(all_datasets)) {
  de_label <- de_label_map[[de_type]]
  
  for (filter_type in names(all_datasets[[de_type]])) {
    dataset_list <- all_datasets[[de_type]][[filter_type]]
    
    for (dataset_name in names(dataset_list)) {
      
      data_dir   <- dataset_list[[dataset_name]]
      output_dir <- file.path(bsj_output_dir, dataset_name)
      
      # Initialize nested list
      if (is.null(limma_results[[dataset_name]])) limma_results[[dataset_name]] <- list()
      if (is.null(limma_results[[dataset_name]][[filter_type]])) limma_results[[dataset_name]][[filter_type]] <- list()
      if (is.null(limma_results[[dataset_name]][[filter_type]][[de_type]])) limma_results[[dataset_name]][[filter_type]][[de_type]] <- list()
      
      message("\n========================================")
      message("Processing: ", de_type, " | ", filter_type, " | ", dataset_name)
      message("========================================")
      
      ## ---- Quantile (once per dataset) ----
      iter_start <- Sys.time()  # Timer for quantile
      message("Running limma quantile | ", dataset_name, " | ", filter_type, " | ", de_type)
      
      # Skip check for quantile
      quantile_file <- sprintf("%s_%s_%s_limma_quantile_summary.csv", dataset_name, filter_type, de_label)
      
      if (file.exists(file.path(output_dir, quantile_file))) {
        message("  Quantile: Skipped, already done.")
        quantile_res <- list(
          metrics = read.csv(file.path(output_dir, quantile_file)),
          de_results = NULL
        )
      } else {
        quantile_res <- run_limma_quantile(
          data_dir = data_dir,
          dataset_name = dataset_name,
          filter_type = filter_type,
          de_label = de_label,
          contrast_formula = contrast_map[[dataset_name]],
          output_dir = output_dir
        )
      }
      
      iter_time <- difftime(Sys.time(), iter_start, units = "secs")
      message("  Quantile completed in ", round(iter_time, 2), " seconds")
      
      limma_results[[dataset_name]][[filter_type]][[de_type]]$quantile <- list(
        result     = if (is.list(quantile_res) && "metrics" %in% names(quantile_res)) quantile_res$metrics else quantile_res,
        de_results = if (is.list(quantile_res) && "de_results" %in% names(quantile_res)) quantile_res$de_results else NULL,
        runtime_seconds = as.numeric(iter_time)
      )
      
      ## ---- Voom (per variant and norm) ----
      limma_results[[dataset_name]][[filter_type]][[de_type]]$voom <- list()
      
      for (variant in variants) {
        for (norm in norm_methods) {
          iter_start <- Sys.time()  # Timer for each voom iteration
          
          message("Running limma voom | ", dataset_name, " | ", filter_type, " | ",
                  de_type, " | ", variant, " | ", norm)
          
          # Skip check for voom
          voom_file <- sprintf("%s_%s_%s_limma_%s_%s_summary.csv", dataset_name, filter_type, de_label, variant, norm)
          
          if (file.exists(file.path(output_dir, voom_file))) {
            message("  Voom ", variant, "/", norm, ": Skipped, already done.")
            voom_res <- list(
              metrics = read.csv(file.path(output_dir, voom_file)),
              de_results = NULL
            )
          } else {
            voom_res <- run_limma_voom(
              data_dir = data_dir,
              dataset_name = dataset_name,
              filter_type = filter_type,
              de_label = de_label,
              norm = norm,
              contrast_formula = contrast_map[[dataset_name]],
              output_dir = output_dir,
              variant = variant
            )
          }
          
          iter_time <- difftime(Sys.time(), iter_start, units = "secs")
          message("  Completed in ", round(iter_time, 2), " seconds")
          
          limma_results[[dataset_name]][[filter_type]][[de_type]]$voom[[variant]][[norm]] <- list(
            result     = if (is.list(voom_res) && "metrics" %in% names(voom_res)) voom_res$metrics else voom_res,
            de_results = if (is.list(voom_res) && "de_results" %in% names(voom_res)) voom_res$de_results else NULL,
            runtime_seconds = as.numeric(iter_time)
          )
        }
      }
    }
  }
}

limma_end_time <- Sys.time()
limma_runtime <- difftime(limma_end_time, limma_start, units = "mins")

# Add overall runtime to results
limma_results[["total_runtime_mins"]] <- as.numeric(limma_runtime)

message(
  "\n========================================\n",
  "LIMMA-VOOM ANALYSIS COMPLETE\n",
  "========================================\n",
  "Total runtime: ", round(limma_runtime, 2), " minutes.\n",
  "========================================"
)

######### EDGER loop ##########
edger_start <- Sys.time()
edger_results <- list()

for (de_type in names(all_datasets)) {
  de_label <- de_label_map[[de_type]]
  
  for (filter_type in names(all_datasets[[de_type]])) {
    dataset_list <- all_datasets[[de_type]][[filter_type]]
    
    for (dataset_name in names(dataset_list)) {
      data_dir <- dataset_list[[dataset_name]]
      output_dir <- file.path(bsj_output_dir, dataset_name)
      
      # Initialize nested list
      if (is.null(edger_results[[dataset_name]])) edger_results[[dataset_name]] <- list()
      if (is.null(edger_results[[dataset_name]][[filter_type]])) edger_results[[dataset_name]][[filter_type]] <- list()
      if (is.null(edger_results[[dataset_name]][[filter_type]][[de_type]])) edger_results[[dataset_name]][[filter_type]][[de_type]] <- list()
      
      file_pattern <- paste0(dataset_name, "_", filter_type, "_counts_", de_label, "_rep")
      
      message("\n========================================")
      message("Processing: ", de_type, " | ", filter_type, " | ", dataset_name)
      message("========================================")
      
      for (norm in norm_methods) {
        iter_start <- Sys.time()
        message("Running EdgeR | ", de_type, " | ", filter_type, " | ", dataset_name, " | ", norm)
        
        # Choose which function to call based on DE type
        edger_res <- if (de_label == "DE0") {
          analyze_edger_DE0(
            data_dir = data_dir,
            dataset_name = dataset_name,
            file_pattern = file_pattern,
            norm_method = norm,
            output_dir = output_dir,
            contrast_formula = contrast_map[[dataset_name]],
            filter = filter_type
          )
        } else {
          analyze_simulation_batch_with_truth(
            data_dir = data_dir,
            dataset_name = dataset_name,
            file_pattern = file_pattern,
            runner_fun = run_edgeR_vdefault_with_truth,   
            norm_method = norm,
            contrast_formula = contrast_map[[dataset_name]],
            threshold_col = "FDR",
            output_dir = output_dir,
            filter_type = filter_type,
            tool = "edgeR"
          )
        }
        
        iter_time <- difftime(Sys.time(), iter_start, units = "secs")
        message("  Completed in ", round(iter_time, 2), " seconds")
        
        # Store result per norm
        edger_results[[dataset_name]][[filter_type]][[de_type]][[norm]] <- list(
          result = if (is.list(edger_res) && "metrics" %in% names(edger_res)) edger_res$metrics else edger_res,
          de_results = if (is.list(edger_res) && "de_results" %in% names(edger_res)) edger_res$de_results else NULL,
          runtime_seconds = as.numeric(iter_time)
        )
        
      }
    }
  }
}

edger_runtime <- difftime(Sys.time(), edger_start, units = "mins")
edger_results[["total_runtime_mins"]] <- as.numeric(edger_runtime)

message(
  "\n========================================\n",
  "EDGER ANALYSIS COMPLETE\n",
  "========================================\n",
  "Total runtime: ", round(edger_runtime, 2), " minutes.\n",
  "========================================"
)

###### DESeq2 loop #####
deseq_start <- Sys.time()
deseq_results <- list()

for (de_type in names(all_datasets)) {
  de_label <- de_label_map[[de_type]]
  
  for (filter_type in names(all_datasets[[de_type]])) {
    dataset_list <- all_datasets[[de_type]][[filter_type]]
    
    for (dataset_name in names(dataset_list)) {
      data_dir <- dataset_list[[dataset_name]]
      output_dir <- file.path(bsj_output_dir, dataset_name)
      
      # Initialize nested list
      if (is.null(deseq_results[[dataset_name]])) deseq_results[[dataset_name]] <- list()
      if (is.null(deseq_results[[dataset_name]][[filter_type]])) deseq_results[[dataset_name]][[filter_type]] <- list()
      if (is.null(deseq_results[[dataset_name]][[filter_type]][[de_type]])) deseq_results[[dataset_name]][[filter_type]][[de_type]] <- list()
      
      file_pattern <- paste0(dataset_name, "_", filter_type, "_counts_", de_label, "_rep")
      file_output_dir <- file.path(output_dir, dataset_name)
      
      message("\n========================================")
      message("Processing: ", de_type, " | ", filter_type, " | ", dataset_name)
      message("========================================")
      
      # Loop over DESeq2 methods
      for (method in c("BetaPrior", "LRT", "WaldTest")) {
        iter_start <- Sys.time()
        
        message("Running DESeq2 | ", method, " | ", de_type, " | ", filter_type, " | ", dataset_name)
        
        output_file <- sprintf("%s_%s_%s_DESeq2_%s_summary.csv", dataset_name, filter_type, de_label, method)
        
        # Skip if already exists
        if (file.exists(file.path(file_output_dir, output_file))) {
          message("Skipped, already done.")
          method_res <- read.csv(file.path(file_output_dir, output_file))
        } else {
          method_res <- switch(method,
                               BetaPrior = if (de_label == "DE0") {
                                 DESeq2_betaprior_DE0(
                                   data_dir = data_dir,
                                   dataset_name = dataset_name,
                                   file_pattern = file_pattern,
                                   output_dir = output_dir,
                                   contrast_formula = contrast_map[[dataset_name]],
                                   filter = filter_type
                                 )
                               } else {
                                 analyze_simulation_batch_with_truth(
                                   data_dir = data_dir,
                                   dataset_name = dataset_name,
                                   file_pattern = file_pattern,
                                   runner_fun = run_DESeq2_BetaPriorT_with_truth, # your DE10p runner
                                   contrast_formula = contrast_map[[dataset_name]],
                                   threshold_col = "padj",
                                   output_dir = output_dir,
                                   filter_type = filter_type,
                                   tool = "DESeq2",
                                   method = method
                                 )
                               },
                               LRT = if (de_label == "DE0") {
                                 analyze_DESeq2_LRT_DE0(
                                   data_dir = data_dir,
                                   dataset_name = dataset_name,
                                   file_pattern = file_pattern,
                                   output_dir = output_dir,
                                   contrast_formula = contrast_map[[dataset_name]],
                                   filter = filter_type
                                 )
                               } else {
                                 analyze_simulation_batch_with_truth(
                                   data_dir = data_dir,
                                   dataset_name = dataset_name,
                                   file_pattern = file_pattern,
                                   runner_fun = run_DESeq2_LRT_with_truth,
                                   contrast_formula = contrast_map[[dataset_name]],
                                   threshold_col = "padj",
                                   output_dir = output_dir,
                                   filter_type = filter_type,
                                   tool = "DESeq2",
                                   method = method
                                 )
                               },
                               WaldTest = if (de_label == "DE0") {
                                 analyze_deseq2_wald_null(
                                   data_dir = data_dir,
                                   dataset_name = dataset_name,
                                   file_pattern = file_pattern,
                                   output_dir = output_dir,
                                   contrast_formula = contrast_map[[dataset_name]],
                                   filter = filter_type
                                 )
                               } else {
                                 analyze_simulation_batch_with_truth(
                                   data_dir = data_dir,
                                   dataset_name = dataset_name,
                                   file_pattern = file_pattern,
                                   runner_fun = run_DESeq2_WaldTest_with_truth,
                                   contrast_formula = contrast_map[[dataset_name]],
                                   threshold_col = "padj",
                                   output_dir = output_dir,
                                   filter_type = filter_type,
                                   tool = "DESeq2",
                                   method = method
                                 )
                               }
          )
        }
        
        iter_time <- difftime(Sys.time(), iter_start, units = "secs")
        message("  Completed in ", round(iter_time, 2), " seconds")
        
        # Store method result
        deseq_results[[dataset_name]][[filter_type]][[de_type]][[method]] <- list(
          result = if (is.list(method_res) && "metrics" %in% names(method_res)) method_res$metrics else method_res,
          de_results = if (is.list(method_res) && "de_results" %in% names(method_res)) method_res$de_results else NULL,
          runtime_seconds = as.numeric(iter_time)
        )
      }
      
      dataset_time <- difftime(Sys.time(), iter_start, units = "mins")
      message("Dataset ", dataset_name, " completed in ", round(dataset_time, 2), " minutes\n")
    }
  }
}

deseq_runtime <- difftime(Sys.time(), deseq_start, units = "mins")
deseq_results[["total_runtime_mins"]] <- as.numeric(deseq_runtime)

message(
  "\n========================================\n",
  "DESEQ2 ANALYSIS COMPLETE\n",
  "========================================\n",
  "Total runtime: ", round(deseq_runtime, 2), " minutes.\n",
  "========================================"
)

###### Run real data ######
all_results <- list()
runtime_log <- list()
plot_out_base <- "/media/meteor/FatDawg/Benchmark_Paper/DE_Results/BSJ-FSJ/Dispersion_plots"

for (dataset_name in names(bsj_paths)) {
  cat("\n========================================\n")
  cat(sprintf("Processing dataset: %s\n", dataset_name))
  cat("========================================\n")
  
  dataset_start <- proc.time()
  
  bsj_path    <- bsj_paths[[dataset_name]]
  fsj_path    <- fsj_paths[[dataset_name]]
  starFC_path <- star_fc_paths[[dataset_name]]
  meta_path   <- metadata_paths[[dataset_name]]
  
  plot_out_dir <- file.path(plot_out_base, dataset_name)
  if (!dir.exists(plot_out_dir)) dir.create(plot_out_dir, recursive = TRUE)
  
  all_results[[dataset_name]] <- list()
  
  log_time <- function(tool, norm, mode, min_cpm, elapsed) {
    runtime_log[[length(runtime_log) + 1]] <<- data.frame(
      dataset     = dataset_name,
      tool        = tool,
      norm        = norm,
      mode        = mode,
      min_cpm     = min_cpm,
      elapsed_sec = round(elapsed, 3),
      stringsAsFactors = FALSE
    )
  }
  
  # ── CiriDE ──────────────────────────────────────────────────────────────────
  for (norm in norm_methods) {
    cat("Running CiriDE-FSJ", norm, "...\n")
    t0 <- proc.time()
    all_results[[dataset_name]][["ciriDE"]][["ciriFSJ"]][[norm]] <- tryCatch({
      run_ciriDE(bsj_path, fsj_path, meta_path, norm, contrast_map[[dataset_name]])
    }, error = function(e) {
      cat(sprintf("  ERROR in CiriDE-FSJ for %s: %s\n", dataset_name, e$message)); NULL
    })
    save_plots(all_results[[dataset_name]][["ciriDE"]][["ciriFSJ"]][[norm]]$plots,
               plot_out_dir, paste0(dataset_name, "_CiriDE_FSJ_", norm))
    log_time("CiriDE", norm, "FSJ", NA, (proc.time() - t0)["elapsed"])
    cat(sprintf("  CiriDE-FSJ %s done in %.2f sec\n", norm, (proc.time() - t0)["elapsed"]))
    
    cat("Running CiriDE-STAR", norm, "...\n")
    t0 <- proc.time()
    all_results[[dataset_name]][["ciriDE"]][["STAR"]][[norm]] <- tryCatch({
      run_ciriDE(bsj_path, starFC_path, meta_path, norm, contrast_map[[dataset_name]])
    }, error = function(e) {
      cat(sprintf("  ERROR in CiriDE-STAR for %s: %s\n", dataset_name, e$message)); NULL
    })
    save_plots(all_results[[dataset_name]][["ciriDE"]][["STAR"]][[norm]]$plots,
               plot_out_dir, paste0(dataset_name, "_CiriDE_STAR_", norm))
    log_time("CiriDE", norm, "STAR", NA, (proc.time() - t0)["elapsed"])
    cat(sprintf("  CiriDE-STAR %s done in %.2f sec\n", norm, (proc.time() - t0)["elapsed"]))
    
    for (min_cpm in c("def", "1", "5")) {
      cat("Running ciriDE-BSJonly", norm, min_cpm, "...\n")
      t0 <- proc.time()
      all_results[[dataset_name]][["ciriDE"]][["BSJonly"]][[norm]][[min_cpm]] <- tryCatch({
        ciriDE_BSJ(bsj_path, meta_path, norm, min_cpm, contrast_map[[dataset_name]])
      }, error = function(e) {
        cat(sprintf("  ERROR in ciriDE-BSJonly for %s: %s\n", dataset_name, e$message)); NULL
      })
      save_plots(all_results[[dataset_name]][["ciriDE"]][["BSJonly"]][[norm]][[min_cpm]]$plots,
                 plot_out_dir, paste0(dataset_name, "_ciriDE_BSJonly_", norm, "_cpm", min_cpm))
      log_time("ciriDE", norm, "BSJonly", min_cpm, (proc.time() - t0)["elapsed"])
      cat(sprintf("  ciriDE-BSJonly %s min_cpm=%s done in %.2f sec\n", norm, min_cpm, (proc.time() - t0)["elapsed"]))
    }
  }
  
  # ── edgeR ───────────────────────────────────────────────────────────────────
  for (norm in norm_methods) {
    cat("Running edgeR-FSJ", norm, "...\n")
    t0 <- proc.time()
    all_results[[dataset_name]][["edgeR"]][["ciriFSJ"]][[norm]] <- tryCatch({
      run_edger(bsj_path, fsj_path, meta_path, norm, contrast_map[[dataset_name]])
    }, error = function(e) {
      cat(sprintf("  ERROR in edgeR-ciriFSJ for %s: %s\n", dataset_name, e$message)); NULL
    })
    save_plots(all_results[[dataset_name]][["edgeR"]][["ciriFSJ"]][[norm]]$plots,
               plot_out_dir, paste0(dataset_name, "_edgeR_FSJ_", norm))
    log_time("edgeR", norm, "FSJ", NA, (proc.time() - t0)["elapsed"])
    cat(sprintf("  edgeR-FSJ %s done in %.2f sec\n", norm, (proc.time() - t0)["elapsed"]))
    
    cat("Running edgeR-STAR", norm, "...\n")
    t0 <- proc.time()
    all_results[[dataset_name]][["edgeR"]][["STAR"]][[norm]] <- tryCatch({
      run_edger(bsj_path, starFC_path, meta_path, norm, contrast_map[[dataset_name]])
    }, error = function(e) {
      cat(sprintf("  ERROR in edgeR-STAR for %s: %s\n", dataset_name, e$message)); NULL
    })
    save_plots(all_results[[dataset_name]][["edgeR"]][["STAR"]][[norm]]$plots,
               plot_out_dir, paste0(dataset_name, "_edgeR_STAR_", norm))
    log_time("edgeR", norm, "STAR", NA, (proc.time() - t0)["elapsed"])
    cat(sprintf("  edgeR-STAR %s done in %.2f sec\n", norm, (proc.time() - t0)["elapsed"]))
    
    for (min_cpm in c("def", "1", "5")) {
      cat("Running edgeR-BSJonly", norm, min_cpm, "...\n")
      t0 <- proc.time()
      all_results[[dataset_name]][["edgeR"]][["BSJonly"]][[norm]][[min_cpm]] <- tryCatch({
        edger_bsj(bsj_path, meta_path, norm, min_cpm, contrast_map[[dataset_name]])
      }, error = function(e) {
        cat(sprintf("  ERROR in edgeR-BSJonly for %s: %s\n", dataset_name, e$message)); NULL
      })
      save_plots(all_results[[dataset_name]][["edgeR"]][["BSJonly"]][[norm]][[min_cpm]]$plots,
                 plot_out_dir, paste0(dataset_name, "_edgeR_BSJonly_", norm, "_cpm", min_cpm))
      log_time("edgeR", norm, "BSJonly", min_cpm, (proc.time() - t0)["elapsed"])
      cat(sprintf("  edgeR-BSJonly %s min_cpm=%s done in %.2f sec\n", norm, min_cpm, (proc.time() - t0)["elapsed"]))
    }
  }
  
  # ── Limma-Voom ──────────────────────────────────────────────────────────────
  for (mode in limma_modes) {
    for (norm in norm_methods) {
      cat("Running LV-FSJ", norm, mode, "...\n")
      t0 <- proc.time()
      all_results[[dataset_name]][["LV"]][["ciriFSJ"]][[mode]][[norm]] <- tryCatch({
        run_limma(bsj_path, fsj_path, meta_path, norm, contrast_map[[dataset_name]], run_mode = mode)
      }, error = function(e) {
        cat(sprintf("  ERROR in LV-ciriFSJ for %s: %s\n", dataset_name, e$message)); NULL
      })
      save_plots(all_results[[dataset_name]][["LV"]][["ciriFSJ"]][[mode]][[norm]]$plots,
                 plot_out_dir, paste0(dataset_name, "_limma_FSJ_", mode, "_", norm))
      log_time("limma", norm, paste0("FSJ-", mode), NA, (proc.time() - t0)["elapsed"])
      cat(sprintf("  LV-FSJ %s %s done in %.2f sec\n", norm, mode, (proc.time() - t0)["elapsed"]))
      
      cat("Running LV-STAR", norm, mode, "...\n")
      t0 <- proc.time()
      all_results[[dataset_name]][["LV"]][["STAR"]][[mode]][[norm]] <- tryCatch({
        run_limma(bsj_path, starFC_path, meta_path, norm, contrast_map[[dataset_name]], run_mode = mode)
      }, error = function(e) {
        cat(sprintf("  ERROR in LV-STAR for %s: %s\n", dataset_name, e$message)); NULL
      })
      save_plots(all_results[[dataset_name]][["LV"]][["STAR"]][[mode]][[norm]]$plots,
                 plot_out_dir, paste0(dataset_name, "_limma_STAR_", mode, "_", norm))
      log_time("limma", norm, paste0("STAR-", mode), NA, (proc.time() - t0)["elapsed"])
      cat(sprintf("  LV-STAR %s %s done in %.2f sec\n", norm, mode, (proc.time() - t0)["elapsed"]))
      
      for (min_cpm in c("def", "1", "5")) {
        cat("Running LV-BSJonly", norm, mode, min_cpm, "...\n")
        t0 <- proc.time()
        all_results[[dataset_name]][["LV"]][["BSJonly"]][[mode]][[norm]][[min_cpm]] <- tryCatch({
          limma_bsj(bsj_path, meta_path, norm, min_cpm, contrast_map[[dataset_name]], run_mode = mode)
        }, error = function(e) {
          cat(sprintf("  ERROR in LV-BSJonly for %s: %s\n", dataset_name, e$message)); NULL
        })
        save_plots(all_results[[dataset_name]][["LV"]][["BSJonly"]][[mode]][[norm]][[min_cpm]]$plots,
                   plot_out_dir, paste0(dataset_name, "_limma_BSJonly_", mode, "_", norm, "_cpm", min_cpm))
        log_time("limma", norm, paste0("BSJonly-", mode), min_cpm, (proc.time() - t0)["elapsed"])
        cat(sprintf("  LV-BSJonly %s %s min_cpm=%s done in %.2f sec\n", norm, mode, min_cpm, (proc.time() - t0)["elapsed"]))
      }
    }
  }
  
  # ── DESeq2 ──────────────────────────────────────────────────────────────────
  for (mode in deseq_modes) {
    cat("Running DESeq2-FSJ", mode, "...\n")
    t0 <- proc.time()
    all_results[[dataset_name]][["DESeq2"]][["ciriFSJ"]][[mode]] <- tryCatch({
      run_deseq(bsj_path, fsj_path, meta_path, "TMM", contrast_map[[dataset_name]], run_mode = mode)
    }, error = function(e) {
      cat(sprintf("  ERROR in DESeq2-FSJ for %s: %s\n", dataset_name, e$message)); NULL
    })
    save_plots(all_results[[dataset_name]][["DESeq2"]][["ciriFSJ"]][[mode]]$plots,
               plot_out_dir, paste0(dataset_name, "_DESeq2_FSJ_", mode))
    log_time("DESeq2", "TMM", paste0("FSJ-", mode), NA, (proc.time() - t0)["elapsed"])
    cat(sprintf("  DESeq2-FSJ %s done in %.2f sec\n", mode, (proc.time() - t0)["elapsed"]))
    
    cat("Running DESeq2-STAR", mode, "...\n")
    t0 <- proc.time()
    all_results[[dataset_name]][["DESeq2"]][["STAR"]][[mode]] <- tryCatch({
      run_deseq(bsj_path, starFC_path, meta_path, "TMM", contrast_map[[dataset_name]], run_mode = mode)
    }, error = function(e) {
      cat(sprintf("  ERROR in DESeq2-STAR for %s: %s\n", dataset_name, e$message)); NULL
    })
    save_plots(all_results[[dataset_name]][["DESeq2"]][["STAR"]][[mode]]$plots,
               plot_out_dir, paste0(dataset_name, "_DESeq2_STAR_", mode))
    log_time("DESeq2", "TMM", paste0("STAR-", mode), NA, (proc.time() - t0)["elapsed"])
    cat(sprintf("  DESeq2-STAR %s done in %.2f sec\n", mode, (proc.time() - t0)["elapsed"]))
    
    cat("Running DESeq2-BSJonly", mode, "...\n")
    t0 <- proc.time()
    all_results[[dataset_name]][["DESeq2"]][["BSJonly"]][[mode]] <- tryCatch({
      deseq_bsj(bsj_path, meta_path, contrast_map[[dataset_name]], run_mode = mode)
    }, error = function(e) {
      cat(sprintf("  ERROR in DESeq2-BSJonly for %s: %s\n", dataset_name, e$message)); NULL
    })
    save_plots(all_results[[dataset_name]][["DESeq2"]][["BSJonly"]][[mode]]$plots,
               plot_out_dir, paste0(dataset_name, "_DESeq2_BSJonly_", mode))
    log_time("DESeq2", "TMM", paste0("BSJonly-", mode), NA, (proc.time() - t0)["elapsed"])
    cat(sprintf("  DESeq2-BSJonly %s done in %.2f sec\n", mode, (proc.time() - t0)["elapsed"]))
  }
  
  dataset_elapsed <- (proc.time() - dataset_start)["elapsed"]
  cat(sprintf("\n>>> Dataset %s completed in %.2f sec (%.2f min)\n",
              dataset_name, dataset_elapsed, dataset_elapsed / 60))
}

# Final runtime table
runtime_df <- dplyr::bind_rows(runtime_log)
print(runtime_df)

# Summary by tool
runtime_summary <- runtime_df %>%
  group_by(tool, mode) %>%
  summarise(
    mean_sec  = round(mean(elapsed_sec), 3),
    total_sec = round(sum(elapsed_sec), 3),
    n_runs    = n(),
    .groups   = "drop"
  ) %>%
  arrange(desc(mean_sec))

print(runtime_summary)

## Final metrics of comaprison: Collect Genes, adj.p-vals, logFC direction, Similarity between identified tools
flat_results  <- flatten_results(all_results)
summary_table <- build_summary_table(flat_results)
