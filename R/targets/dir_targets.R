dir_targets = list(
  
  # location of sattag files
  tar_target(sattag_dir, file.path('data', 'sattag')),
  
  # location to store extracted and imputed depth bins
  tar_target(
    name = imputed_out_dir, 
    command = {
      f = file.path('output', 'imputed_bins')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # location to store extracted and imputed depth bins
  tar_target(
    name = posterior_stages_dir, 
    command = {
      f = file.path('output', 'posterior_imputations')
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
  
  # location for storing sattag dive id labels
  tar_target(
    name = dive_label_plot_dir, 
    command = {
      f = file.path('output', 'tag_labels')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # location for storing dive id labels for imputed dives
  tar_target(
    name = imputed_dive_label_plot_dir, 
    command = {
      f = file.path('output', 'tag_labels_imputed')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # location for storing dive segmentation diagnostics
  tar_target(
    name = label_diagnostic_plot_dir, 
    command = {
      f = file.path(dive_label_plot_dir, 'diagnostics')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # location for storing mcmc samples
  tar_target(
    name = mcmc_sample_dir, 
    command = {
      f = file.path('output', 'mcmc')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  )
  
)
