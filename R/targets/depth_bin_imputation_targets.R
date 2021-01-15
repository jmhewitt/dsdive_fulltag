depth_bin_imputation_targets = list(
  
  # extract and impute depth bins from raw data
  tar_target(
    name = imputed_bins_summary,
    command = impute_bins(
      message_files = message_files[!grepl('069', tar_read(depth_files))],
      depth_files = depth_files[!grepl('069', tar_read(depth_files))],
      imputed_out_dir = imputed_out_dir,
      imputed_plot_dir = imputed_plot_dir
    )
  ),
  
  # use the deepest set of depth bins as the template depth bins
  tar_target(
    name = template_bins, 
    command = read.csv(file = as.character(
      imputed_bins_summary %>%
        slice(which.max(max.depth)) %>%
        dplyr::select(file) %>%
        unlist()
    ))
  )
  
)
  