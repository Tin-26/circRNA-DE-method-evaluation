## Samples for simulation ##
datasets <- list(
  `BC` = list("/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/BC_filtered/automated_filtered_BC.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/BC_filtered/min1_filtered_BC.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/BC_filtered/min5_filtered_BC.csv"),
  `EBC1` = list("/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/EBC1_filtered/automatic_filtered_EBC.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/EBC1_filtered/min1_filtered_EBC.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/EBC1_filtered/min5_filtered_EBC.csv"),
  `EBC2` = list("/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/EBC2_filtered/automatic_filtered_EBC2.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/EBC2_filtered/min1_filtered_EBC2.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/EBC2_filtered/min5_filtered_EBC2.csv"),
  `HCC-PBMC` = list("/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/HCC_PBMC_filtered/automated_filtered_HCC_PBMC.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/HCC_PBMC_filtered/min1_filtered_HCC_PBMC.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/HCC_PBMC_filtered/min5_filtered_HCC_PBMC.csv"),
  `HCC-tissue` = list("/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/HCC_Tissue_filtered/automated_filtered_HCC_Tissue.csv","/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/HCC_Tissue_filtered/min1_filtered_HCC_Tissue.csv", "/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/filtered_datasets/HCC_Tissue_filtered/min5_filtered_HCC_Tissue.csv")
)

metadata <- list(
  `BC` = "/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/metadata/PRJNA553624_metadata.csv",
  `EBC1` = "/media/meteor/FatDawg/Benchmark_Paper/Own_dataset/2026-02-19/Metadata_for_HC.csv",
  `EBC2` = "/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/metadata/2025-11-25_EBC2.csv",
  `HCC-PBMC` = "/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/metadata/PBMC_metadata_cleaned.csv",
  `HCC-tissue` = "/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/HC_filtered_final_datasets/metadata/PRJNA716508_metadata_cleaned.csv"
)

hc_filtered <- function(dgelist_filtered,
                        metadata,
                        nclust = 8,
                        Run = "Run",
                        Sample = "Source",
                        filter, dataset) {

  stopifnot(inherits(dgelist_filtered, "DGEList"))
  stopifnot(is.data.frame(metadata) || data.table::is.data.table(metadata))

  # ---- 1) Library sizes from DGEList ----
  lib_sizes <- dgelist_filtered$samples$lib.size
  sample_ids <- rownames(dgelist_filtered$samples)
  if (is.null(sample_ids)) {
    stop("dgelist_filtered$samples must have rownames corresponding to sample IDs (e.g., Run).")
  }
  names(lib_sizes) <- sample_ids

  # ---- 2) Match metadata to DGEList sample order ----
  if (!Run %in% colnames(metadata)) stop(sprintf("Column '%s' not found in metadata.", Run))
  if (!Sample %in% colnames(metadata)) stop(sprintf("Column '%s' not found in metadata.", Sample))

  matched_metadata <- metadata[match(sample_ids, metadata[[Run]]), , drop = FALSE]
  if (anyNA(matched_metadata[[Run]])) {
    missing <- sample_ids[is.na(matched_metadata[[Run]])]
    stop("These DGEList samples were not found in metadata[[Run]]: ",
         paste(missing, collapse = ", "))
  }

  # sanity check
  stopifnot(all(matched_metadata[[Run]] == sample_ids))

  # ---- 3) Build base table ----
  lib_size_dt <- data.table::data.table(
    SampleID = matched_metadata[[Run]],
    Condition = matched_metadata[[Sample]],
    lib_size = as.numeric(lib_sizes)
  )

  # ---- 4) Cluster WITHIN each condition ----
  # If a condition has fewer samples than nclust, reduce k for that condition.
  clust_dt <- lib_size_dt[, {
    x <- lib_size
    k <- min(nclust, length(x))
    hc <- stats::hclust(stats::dist(x, method = "euclidean"), method = "complete")
    .(SampleID = SampleID,
      lib_size = lib_size,
      Condition.x = Condition,
      LibGroup = as.integer(stats::cutree(hc, k = k)))
  }, by = Condition]


  print(clust_dt)

  # ---- 5) Plot (faceted by condition) ----
  plot <- ggplot2::ggplot(
    clust_dt,
    ggplot2::aes(x = factor(LibGroup), y = lib_size)
  ) +
    ggplot2::geom_boxplot(varwidth = TRUE, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, shape = 1) +
    ggplot2::facet_wrap(~Condition.x, scales = "free_x") +
    ggplot2::labs(
      title ="Library Size Clustering (within each condition)",
      subtitle = paste("Dataset: ", dataset, " – filter:", filter),
      x = "Library Size Cluster",
      y = "Library Size (BSJ counts)"
    ) +
    ggplot2::theme_minimal()

  # ---- 6) Summary table ----
  summary_dt <- clust_dt[, .(
    Nsamples = .N,
    MinLibSize = min(lib_size),
    MedianLibSize = median(lib_size),
    MaxLibSize = max(lib_size)
  ), by = .(Condition.x, LibGroup)][order(Condition.x, LibGroup)]
  
  #### New tablef for IQR
  lib_size_sample_summary <- data.frame(libsizes = dgelist_filtered$samples$lib.size,
                                        sample_ids = rownames(dgelist_filtered$samples))

  return(list(
    matched_metadata = matched_metadata,
    sample_libsize_clust_dt = clust_dt,
    plot = plot,
    hc_summary_dt = summary_dt,
    lib_size_sample_summary = lib_size_sample_summary
    ))
}

filter_names <- c("autofilter", "min1", "min5")

datasets <- lapply(datasets, function(x) {
  stopifnot(length(x) == length(filter_names))
  setNames(x, filter_names)
})


all_results <- lapply(names(datasets), function(ds) {
  
  message("Processing dataset: ", ds)
 
  meta <- read.csv(
    metadata[[ds]],
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  ds_filters <- datasets[[ds]]
  
  filter_results <- lapply(names(ds_filters), function(filt) {
    
    message("  Filter: ", filt)

    mtx <- read.csv(
      ds_filters[[filt]],
      sep = ',',
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )
    
    stopifnot(all(colnames(mtx) %in% meta$Sample_name))
    meta <- meta[match(colnames(mtx), meta$Sample_name), ]
    stopifnot(all(meta$Sample_name == colnames(mtx)))
    
    dge <- DGEList(counts = mtx, group=meta$Source)

    hc_res <- hc_filtered(
      dgelist_filtered = dge,
      metadata         = meta,
      nclust           = 8,
      Run              = "Sample_name",
      Sample           = "Source",
      filter           = filt,
      dataset          = ds
    )
    
    list(
      dge      = dge,
      hc       = hc_res,
      n_genes  = nrow(mtx),
      lib_size = dge$samples$lib.size
    )
  })
  
  setNames(filter_results, names(ds_filters))
})

names(all_results) <- names(datasets)

View(all_results)

drop_samples <- list(
  BC = c("SRR11600338"),
  EBC1 = c(),
  EBC2 = c("Healthy7_Run2"),
  `HCC-PBMC` = c("C1", "C3", "C4"),
  `HCC-tissue` = c("SRR14027947", "SRR14027941")
)

filter_and_save_dataset <- function(res, ds_name, outdir) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Determine final sample set per dataset ----
  final_samples_list <- lapply(names(res), function(filt) {
    dge <- res[[filt]]$dge
    
    setdiff(
      colnames(dge$counts),
      drop_samples[[ds_name]]
    )
  })
  
  # Intersection of samples across all filters
  final_samples <- Reduce(intersect, final_samples_list)
  
  # ---- Save metadata once ----
  meta <- res[[1]]$hc$matched_metadata
  meta_filt <- meta[meta$Sample_name %in% final_samples, ]
  write.csv(
    meta_filt,
    file = file.path(outdir, paste0(ds_name, "_metadata.csv")),
    row.names = FALSE
  )
  
  # ---- Filter DGE matrices per filter and save ----
  filtered <- lapply(names(res), function(filt) {
    dge <- res[[filt]]$dge
    dge_filt <- dge[, final_samples, keep.lib.sizes = FALSE]
    
    write.csv(
      dge_filt$counts,
      file = file.path(outdir, paste0(ds_name, "_", filt, "_counts.csv"))
    )
    
    list(
      dge = dge_filt
    )
  })
  
  names(filtered) <- names(res)
  return(filtered)
}

filtered_results <- lapply(names(all_results), function(ds) {
  message("Filtering and saving dataset: ", ds)
  filter_and_save_dataset(
    res     = all_results[[ds]],
    ds_name = ds,
    outdir  = file.path("/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/for_simulation", ds)
  )
})

names(filtered_results) <- names(all_results)
