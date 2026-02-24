library(SPsimSeq)
library(edgeR)
library(SingleCellExperiment)
library(SimSeq)

set.seed(42) # IMPORTANT

simulate_dataset <- function(count_matrix_path,
                             metadata_path,
                             output_path,
                             DE_fraction,
                             dataset_name,
                             filter_name) {
  
  # --- Read data ---
  raw_counts <- read.csv(count_matrix_path, sep = ",", header = TRUE, row.names = 1)
  study_design <- read.csv(metadata_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
  
  stopifnot(all(colnames(raw_counts) %in% study_design$Sample_name))
  study_design <- study_design[match(colnames(raw_counts), study_design$Sample_name), ]
  
  # --- EdgeR DGEList ---
  dge <- DGEList(counts = raw_counts, group = study_design$Source)
  filtered_counts <- dge$counts
  
  # --- Determine tot.samples and n.sim based on dataset size ---
  n_orig <- ncol(filtered_counts)
  if (n_orig <= 12) {
    tot.samples <- 16    # small datasets get 16 samples (8 per group)
    n.sim <- 30
  } else {
    tot.samples <- 2 * floor(n_orig / 2) # size for bigger datasets, must be able to be divided by 2
    n.sim <- 50
  }
  
  # --- Skip if all outputs already exist ---
  expected_files <- unlist(
    lapply(seq_len(n.sim), function(i) {
      prefix <- paste(
        dataset_name,
        filter_name,
        paste0("DE", DE_fraction),
        paste0("rep", i),
        sep = "_"
      )
      c(
        file.path(output_path, paste0(prefix, "_counts.csv")),
        file.path(output_path, paste0(prefix, "_TP.csv")),
        file.path(output_path, paste0(prefix, "_design.csv"))
      )
    })
  )
  
  if (all(file.exists(expected_files))) {
    message("  -> Outputs already exist, skipping simulation.")
    return(invisible(NULL))
  }
  # --- Prepare group factor ---
  grp <- factor(study_design$Source)
  stopifnot(nlevels(grp) == 2)
  group.config <- rep(0.5, 2)
  
  message("Simulating ", n.sim, " datasets for ", dataset_name, " (", filter_name,
          "), DE_fraction = ", DE_fraction, ", total_samples = ", tot.samples)
  
  # --- Run SPsimSeq ---
  sim_data <- SPsimSeq(
    n.sim = n.sim,
    s.data = filtered_counts,
    group = grp,
    group.config = group.config,
    n.genes = min(1000, nrow(filtered_counts)),
    tot.samples = tot.samples,
    pDE = DE_fraction,
    lfc.thrld = 0.5,
    t.thrld = 2.5,
    genewiseCor = TRUE,
    log.CPM.transform = TRUE,
    variable.lib.size = FALSE,
    result.format = "list",
    #return.details = TRUE,
    verbose = TRUE
  )
  sim_data
  cat("Generated", length(sim_data), "datasets.\n")

  # --- Save all simulated datasets ---
  for (i in seq_along(sim_data)) {
    sim <- sim_data[[i]]
    prefix <- paste(dataset_name, filter_name, paste0("DE", DE_fraction), paste0("rep", i), sep = "_")
    
    message("  Saving replicate ", i, "...")
    counts_file <- file.path(output_path, paste0(prefix, "_counts.csv"))
    tp_file     <- file.path(output_path, paste0(prefix, "_TP.csv"))
    design_file <- file.path(output_path, paste0(prefix, "_design.csv"))
    
    if (file.exists(counts_file) &&
        file.exists(tp_file) &&
        file.exists(design_file)) {
      message("  Replicate ", i, " already exists, skipping.")
      next
    }
    write.csv(sim$counts, counts_file)
    write.csv(sim$rowData, tp_file)
    write.csv(sim$colData, design_file, row.names = FALSE)
  }
  
  message("Done with ", dataset_name, " (", filter_name, "), DE_fraction = ", DE_fraction, "\n")
}

# --- Run simulation for all datasets and filters ---
base_outdir <- "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only"
input_folder <- "/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/for_simulation/BSJ_only"
DE_settings <- c(0, 0.1)

datasets <- list.dirs(input_folder, recursive = FALSE)

for (ds_path in datasets) {
  ds_name <- basename(ds_path)
  message("Dataset: ", ds_name)
  
  metadata_path <- file.path(ds_path, paste0(ds_name, "_metadata.csv"))
  
  count_files <- list.files(ds_path, pattern = "_counts\\.csv$", full.names = TRUE)
  
  for (cnt in count_files) {
    filter_name <- sub("\\.csv$", "", sub(paste0(ds_name, "_|_counts"), "", basename(cnt)))
    message("  Filter: ", filter_name)
    
    for (DEf in DE_settings) {
      outdir <- file.path(base_outdir, ds_name, filter_name, paste0("DE_", DEf))
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      simulate_dataset(
        count_matrix_path = cnt,
        metadata_path = metadata_path,
        output_path = outdir,
        DE_fraction = DEf,
        dataset_name = ds_name,
        filter_name = filter_name
      )
    }
  }
}

