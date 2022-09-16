depth_bin_imputation_targets = list(
  
  # location to store extracted and imputed depth bins
  tar_target(
    name = imputed_out_dir, 
    command = {
      f = file.path('output', 'imputed_bins')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # location to store depth bin imputation diagnostics
  tar_target(
    name = imputed_plot_dir, 
    command = {
      f = file.path(imputed_out_dir, 'plots')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # TODO: revise processing in impute_bins b/c the assumptions made do not hold
  #  to larger groups of tags, esp. wrt. how the MaxDepth is used to create 16 
  #  bins.  There seems to be more nuance than just this alone.
  # 
  # # extract and impute depth bins from raw data
  # tar_target(
  #   name = imputed_bins_summary,
  #   command = impute_bins(
  #     message_files = message_files[!grepl('069', tar_read(depth_files))],
  #     depth_files = depth_files[!grepl('069', tar_read(depth_files))],
  #     imputed_out_dir = imputed_out_dir,
  #     imputed_plot_dir = imputed_plot_dir
  #   )
  # ),
  # 
  # # use the deepest set of depth bins as the template depth bins
  # tar_target(
  #   name = template_bins, 
  #   command = read.csv(file = as.character(
  #     imputed_bins_summary %>%
  #       slice(which.max(max.depth)) %>%
  #       dplyr::select(file) %>%
  #       unlist()
  #   ))
  # )
  
  # use the deepest set of depth bins as the template depth bins
  tar_target(
    name = template_bins,
    command = read.csv(
      file = file.path('output', 'imputed_bins', 'bin2352.csv')
    )
  )
  
)
