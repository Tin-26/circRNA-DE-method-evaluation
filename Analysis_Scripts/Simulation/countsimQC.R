library(countsimQC)

base_outdir      <- "/media/meteor/FatDawg/Benchmark_Paper/countsimQC_reports/BSJ_only"
original_inputs  <- "/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/for_simulation/BSJ_only"
simulated_inputs <- "/media/meteor/FatDawg/Benchmark_Paper/Simulated_datasets/BSJ_only"

datasets <- list.dirs(original_inputs, recursive = FALSE, full.names = FALSE)

for (dataset in datasets) {
  
  message("Processing dataset: ", dataset)
  
  orig_dir <- file.path(original_inputs, dataset)
  sim_base <- file.path(simulated_inputs, dataset)
  
  if (!dir.exists(sim_base)) next
  
  ## filter folders: autofilter_counts, min1_counts, ...
  filter_dirs <- list.dirs(sim_base, recursive = FALSE, full.names = FALSE)
  
  for (filter_dir in filter_dirs) {
    print(filter_dir)

    filter_name <- sub("_counts$", "", filter_dir)
    sim_filter_dir <- file.path(sim_base, filter_dir)
    
    ## DE folders: DE_0, DE_0.1, ...
    de_dirs <- list.dirs(sim_filter_dir, recursive = FALSE, full.names = FALSE)
    
    for (de_dir in de_dirs) {
      
      message("  Filter: ", filter_name, " | ", de_dir)
      
      out_dir <- file.path(
        base_outdir,
        dataset,
        filter_name,
        de_dir
      )
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      ## ----------------------------
      ## Original counts
      ## ----------------------------
      orig_file <- file.path(
        orig_dir,
        paste0(dataset, "_", filter_name, "_counts.csv")
      )
      
      if (!file.exists(orig_file)) {
        warning("Missing original file: ", orig_file)
        next
      }
      
      original_matrix <- as.matrix(
        read.csv(orig_file, row.names = 1, check.names = FALSE)
      )
      
      ## ----------------------------
      ## Simulated replicates
      ## ----------------------------
      sim_dir <- file.path(sim_filter_dir, de_dir)
      
      sim_files <- list.files(
        sim_dir,
        pattern = "_counts\\.csv$",
        full.names = TRUE
      )
      
      if (length(sim_files) == 0) {
        warning("No simulated files in ", sim_dir)
        next
      }
      
      count_list <- lapply(sim_files, function(f) {
        as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
      })
      
      names(count_list) <- sub("_counts\\.csv$", "", basename(sim_files))
      
      count_list[["original"]] <- original_matrix

      ## ----------------------------
      ## Run countsimQC
      ## ----------------------------
      countsimQCReport(
        ddsList     = count_list,
        outputFile = "countsimQC_report.html",
        outputDir  = out_dir,
        savePlots  = TRUE,
        description = paste(
          "Dataset:", dataset,
          "| Filter:", filter_name,
          "|", de_dir
        ),
        dpi = 300
      )
    }
  }
}
