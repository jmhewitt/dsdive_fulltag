report_targets = list(

  # list of directories containing mcmc output
  tarchetypes::tar_files(sample_dirs, dir(mcmc_sample_dir, full.names = TRUE)),
  
  # convergence and sampler diagnostics for model parameters
  tarchetypes::tar_render_rep(
    name = parameter_diagnostics,
    path = 'reports/parameter_diagnostics/parameter_diagnostics.Rmd',
    params = data.frame(
      tgt_dir = basename(sample_dirs),
      output_file = paste(
        'parameter_diagnostics_',
        basename(sample_dirs),
        '.pdf',
        sep = ''
      )
    ),
    deployment = 'worker'
  ),
  
  # posterior plots for stages
  tar_target(
    name = posterior_imputation_plots,
    command = imputation_plots(
      mcmc_sample_dir = sample_dirs,
      output_dir = file.path(posterior_stages_dir, basename(sample_dirs)),
      template_bins = template_bins,
      tag_info = tag_info
    ), 
    pattern = map(sample_dirs),
    deployment = 'worker', 
  )
  
)
