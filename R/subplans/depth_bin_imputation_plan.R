depth_bin_imputation_plan = drake_plan(
  
  # location for storing extracted and imputed depth bins
  imputed_out_dir = file_out(!!file.path('data', 'imputed_bins')),
  
  # location for storing diagnostic plots of imputed depth bins
  imputed_plot_dir = file_out(!!file.path('plots', 'imputed_bins')),
  
  # extract and impute depth bins from raw data
  imputed_bins_summary = impute_bins(
    message_files = message_files,
    depth_files = depth_files,
    imputed_out_dir = imputed_out_dir,
    imputed_plot_dir = imputed_plot_dir
  ),
  
  # use the deepest set of depth bins as the template depth bins
  template_bins = read.csv(file = as.character(
    imputed_bins_summary %>% 
      slice(which.max(max.depth)) %>% 
      dplyr::select(file) %>% 
      unlist()
  ))
  
)