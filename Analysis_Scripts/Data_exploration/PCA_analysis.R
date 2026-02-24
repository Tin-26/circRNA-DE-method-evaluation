# Visualize PCAs using the BSJ, FSJ, and FSJ and FSJ count matrices

# Function 1: Generates the PCA plots for the BSJ and FSJ feature counts

run_bsj_pca <- function(
    metadata_path,
    bsj_path,
    circ_id_col = "circ_id",
    gene_col = "geneName",
    metadata_sample_col = "Sample_ID",
    metadata_group_col = "Disease",
    title = "BSJ-Only PCA Plot"
) {
  
  # ---- Load libraries ----
  library(readxl)
  library(readr)
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(tibble)
  
  # ---- Load the study metadata ----
  message("Reading metadata...")
  metadata <- if (grepl("\\.xlsx$", metadata_path)) {
    read_excel(metadata_path)
  } else {
    read_csv(metadata_path)
  }
  
  # ---- Load the BSJ or FSJ matrix ----
  message("Reading the BSJ or FSJ matrix...")
  bsj_mtx <- read_table(bsj_path, show_col_types = FALSE)
  
  # ---- Remove annotation columns ----
  # not needed at this stage
  bsj_counts_mtx <- bsj_mtx %>%
    dplyr::select(-any_of(c(circ_id_col, gene_col)))
  
  # ---- Metadata matching ----
  # very important to not mix the samples and their conditions
  metadata_matched <- metadata[
    match(colnames(bsj_counts_mtx), metadata[[metadata_sample_col]]),
  ]
  
  checkpoint_1 <- all(colnames(bsj_counts_mtx) == metadata_matched[[metadata_sample_col]])
  if (!checkpoint_1) {
    message("Column names present in the BSJ matrix:")
    print(colnames(bsj_counts_mtx))
    message("Sample IDs in metadata:")
    print(metadata[[metadata_sample_col]])
    stop("Checkpoint 1 failed: BSJ columns do NOT match metadata sample IDs.")
  } else {
    message("Checkpoint 1 passed: BSJ column order matches metadata.")
  }
  
  # ---- Define group ----
  group <- factor(metadata_matched[[metadata_group_col]])
  
  # ---- Build DGEList ----
  dgelist <- DGEList(counts = bsj_counts_mtx, group = group)
  
  # ---- Checkpoint 2: DGEList rownames vs metadata ----
  samples <- dgelist$samples
  
  checkpoint_2 <- all(rownames(samples) == metadata_matched[[metadata_sample_col]]) &&
    all(samples$group == metadata_matched[[metadata_group_col]])
  
  if (!checkpoint_2) {
    stop("Checkpoint 2 failed: DGEList samples do NOT match metadata.")
  } else {
    message("Checkpoint 2 passed: DGEList rownames & group match metadata.")
  }
  
  # ---- Filter low expression circRNAs ----
  keep <- filterByExpr(dgelist, group = group)
  
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # ---- Checkpoint 3: filtered DGEList ----
  samples_filtered <- dgelist_filtered$samples
  
  checkpoint_3 <- all(rownames(samples_filtered) == metadata_matched[[metadata_sample_col]]) &&
    all(samples_filtered$group == metadata_matched[[metadata_group_col]])
  
  if (!checkpoint_3) {
    stop("Checkpoint 3 failed: Filtered DGEList does NOT match metadata.")
  } else {
    message("Checkpoint 3 passed: Filtered DGEList matches metadata.")
  }
  
  # ---- Normalize ----
  dgelist_norm <- calcNormFactors(dgelist_filtered, method = "TMM")
  
  # ---- Generate the PCA ----
  log2.cpm <- cpm(dgelist_norm, log = TRUE)
  pca.res <- prcomp(t(log2.cpm), scale. = FALSE)
  
  pca.res.df <- as_tibble(pca.res$x)
  pca.res.df$Group <- group
  
  # reorder metadata to PCA order
  metadata_ordered <- metadata_matched[
    match(rownames(pca.res$x), metadata_matched[[metadata_sample_col]]),
  ]
  
  # ---- Checkpoint 4: PCA rownames vs metadata ----
  checkpoint_4 <- all(rownames(pca.res$x) == metadata_ordered[[metadata_sample_col]])
  
  if (!checkpoint_4) {
    stop("Checkpoint 4 failed: PCA rownames do NOT match metadata order.")
  } else {
    message("Checkpoint 4 passed: PCA sample order matches metadata.")
  }
  
  pca.res.df$Sample <- metadata_ordered[[metadata_sample_col]]
  
  # ---- Variance ----
  pc.var <- pca.res$sdev^2
  pc.per <- round(pc.var / sum(pc.var) * 100, 1)
  
  # ---- PCA plot ----
  pca_plot <- ggplot(pca.res.df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 6, shape = 16) +
    xlab(paste0("PC1 (", pc.per[1], "%)")) +
    ylab(paste0("PC2 (", pc.per[2], "%)")) +
    labs(title = title) +
    theme_classic()
  
  # ---- Return ----
  return(list(
    pca_plot = pca_plot,
    pca_df = pca.res.df,
    dgelist = dgelist_norm
  ))
}

# Run function requires:
# 1. Path to the study metadata
# 2. Path to the BSJ feature counts matrix
# 3. Title of the column that contains the circRNA IDs (can be left empty with "" if your data does not contain this column)
# 4. Title of the column that contains the host gene names (can be left empty with "" if your data does not contain this column)
# 5. Title of the column in your metadata file that matches the colnames of the BSJ and FSJ matrix
# 6. Title of the column in your metadata file that specifies the conditions being tested.
# 7. Title of the plot (optional, can be left empty with "" if a title is not needed).

res <- run_bsj_pca(
  metadata_path = "path/to/metadata_csv/metadata.csv",
  bsj_path = "path/to/bsj_mtx/bsj_mtx.txt",
  circ_id_col = "circ_id",
  gene_col = "GeneName",
  metadata_sample_col = "Sample Name",
  metadata_group_col = "Group",
  title = "Test Title"
)

# The output is a list which contains the following information:
# 1. The dataframe used to generate the PCA
# 2. The generated PCA plot
# 3. BSJ/FSJ DGEList
# 4. BSJ/FSJ filtered DGEList
# 5. BSJ/FSJ filtered and normalized DGEList

# Access the generated PCA plot
res$pca_plot


############## Function 2: Generates the PCA plots for the combined FSJ and BSJ ############## 

run_fsj_bsj_pca <- function(
    metadata_path,
    bsj_path,
    fsj_path,
    circ_id_col = "circRNA_ID",
    gene_col = "gene_id",
    metadata_sample_col = "Run",
    metadata_group_col = "Source",
    title = ""
) {
  
  # ---- Load libraries ----
  library(readxl)
  library(readr)
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(tibble)
  
  # ---- Load metadata ----
  message("Reading metadata...")
  metadata <- if (grepl("\\.xlsx$", metadata_path)) {
    read_excel(metadata_path)
  } else {
    read_csv(metadata_path)
  }
  
  # ---- Load BSJ matrix ----
  message("Reading the BSJ matrix")
  bsj_mtx <- read.table(bsj_path, sep = "\t", header = TRUE)
  
  # ---- Load FSJ matrix ----
  message("Reading the FSJ matrix")
  fsj_mtx <- read.table(fsj_path, sep = "\t", header = TRUE)
  
  # ---- Remove annotation columns ----
  remove_cols <- c(circ_id_col, gene_col)
  
  remove_cols_bsj <- remove_cols[remove_cols %in% colnames(bsj_mtx)]
  remove_cols_fsj <- remove_cols[remove_cols %in% colnames(fsj_mtx)]
  
  bsj_counts_mtx <- bsj_mtx %>%
    dplyr::select(-any_of(remove_cols_bsj))
  
  fsj_counts_mtx <- fsj_mtx %>%
    dplyr::select(-any_of(remove_cols_fsj))
  
  # ---- Match metadata to FSJ ----
  metadata_matched_fsj <- metadata[
    match(colnames(fsj_counts_mtx), metadata[[metadata_sample_col]]),
  ]
  
  stopifnot(all(colnames(fsj_counts_mtx) == metadata_matched_fsj[[metadata_sample_col]]))
  
  # ---- Match metadata to BSJ ----
  metadata_matched_bsj <- metadata[
    match(colnames(bsj_counts_mtx), metadata[[metadata_sample_col]]),
  ]
  
  stopifnot(all(colnames(bsj_counts_mtx) == metadata_matched_bsj[[metadata_sample_col]]))
  
  # ---- Define group ----
  group <- factor(metadata_matched_fsj[[metadata_group_col]])
  
  # ---- Build FSJ DGEList ----
  fsj_dgelist <- DGEList(counts = fsj_counts_mtx, group = group)
  
  # ---- Sanity check 1 ----
  samples <- fsj_dgelist$samples
  
  condition_match <- all(rownames(samples) == metadata_matched_fsj[[metadata_sample_col]]) &&
    all(samples$group == metadata_matched_fsj[[metadata_group_col]])
  
  if (condition_match) {
    message("DGEList rownames AND metadata group both match - continue the analysis.")
  } else {
    stop("Mismatch detected: rownames or group do NOT match metadata.")
  }
  
  # ---- Filter low expression ----
  fsj_keep <- filterByExpr(fsj_dgelist, group = group)
  
  fsj_dgelist_filtered <- fsj_dgelist[fsj_keep, , keep.lib.sizes = FALSE]
  
  # ---- Sanity check 2 ----
  samples_filtered <- fsj_dgelist_filtered$samples
  
  condition_match_filtered <- all(rownames(samples_filtered) == metadata_matched_fsj[[metadata_sample_col]]) &&
    all(samples$group == metadata_matched_fsj[[metadata_group_col]])
  
  if (condition_match_filtered) {
    message("Filtered DGEList rownames AND metadata group both match - continue the analysis.")
  } else {
    stop("Mismatch detected: rownames or group do NOT match metadata.")
  }
  
  # ---- Normalize FSJ ----
  fsj_dgelist_filtered_norm <- calcNormFactors(fsj_dgelist_filtered, method = "TMM")
  
  # ---- Apply FSJ normalization to BSJ ----
  message("Applying FSJ normalization to BSJ...")
  
  circRNA_DGE <- DGEList(
    counts = bsj_counts_mtx,
    group = group,
    lib.size = fsj_dgelist_filtered_norm$samples$lib.size,
    norm.factors = fsj_dgelist_filtered_norm$samples$norm.factors
  )
  
  # ---- Sanity check 3 ----
  circRNA_DGE_samples <- circRNA_DGE$samples
  
  condition_match_circRNA_DGE <- all(rownames(circRNA_DGE_samples) == metadata_matched_fsj[[metadata_sample_col]]) &&
    all(circRNA_DGE_samples$group == metadata_matched_fsj[[metadata_group_col]])
  
  if (condition_match_circRNA_DGE) {
    message("circRNA DGEList rownames AND metadata group both match - continue the analysis.")
  } else {
    stop("Mismatch detected: circRNA DGEList rownames or group do NOT match metadata.")
  }
  
  # ---- PCA ----
  log2.cpm <- cpm(circRNA_DGE, log = TRUE)
  
  pca.res <- prcomp(t(log2.cpm), scale. = FALSE)
  
  pc.var <- pca.res$sdev^2
  pc.per <- round(pc.var / sum(pc.var) * 100, 1)
  
  # ---- PCA dataframe ----
  pca.res.df <- as_tibble(pca.res$x)
  pca.res.df$Group <- group
  pca.res.df$Sample <- metadata_matched_fsj[[metadata_sample_col]]
  
  metadata_matched_fsj_ordered <- metadata_matched_fsj[
    match(rownames(pca.res$x), metadata_matched_fsj[[metadata_sample_col]]),
  ]
  
  pca.res.df$Sample <- metadata_matched_fsj_ordered[[metadata_sample_col]]
  
  # ---- Sanity check 4 ----
  pca.res.x <- pca.res$x
  
  condition_match_pca.res.x <- all(rownames(pca.res.x) == pca.res.df$Sample)
  
  if (condition_match_pca.res.x) {
    message("pca.res.x rownames AND pca.res.df sample order match - continue the analysis.")
  } else {
    stop("Mismatch detected: pca.res.x rownames and pca.res.df sample order DO NOT match.")
  }
  
  # ---- PCA plot ----
  pca_plot <- ggplot(pca.res.df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 6, shape = 16) +
    xlab(paste0("PC1 (", pc.per[1], "%)")) +
    ylab(paste0("PC2 (", pc.per[2], "%)")) +
    labs(title = title) +
    theme_classic()
  
  
  return(list(
    pca_plot = pca_plot,
    pca_df = pca.res.df,
    fsj_dge = fsj_dgelist,
    fsj_dgelist_filtered = fsj_dgelist_filtered,
    fsj_dgelist_filtered_norm = fsj_dgelist_filtered_norm
  ))
}

res <- run_fsj_bsj_pca(
  metadata_path = "/path/to/metadata_files.csv",
  bsj_path = "/path/to/bsj_feature_counts_file/bsj_feature_counts.txt",
  fsj_path = "/path/to/fsj_feature_counts_file/fsj_feature_counts.txt",
  circ_id_col = "circ_id",
  gene_col = "GeneName",
  metadata_sample_col = "Sample Name",
  metadata_group_col = "Group",
  title = "Test Title"
)
