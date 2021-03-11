validation_targets = list(
  
  tar_target(
    name = nim_pkg_val, 
    command = flatten_tags(
      tag_list = imputed_dive, 
      transition_matrices = transition_matrices, 
      template_bins = template_bins, 
      movement_types = movement_types, 
      pi_discretization = parameter_discretization$pi, 
      lambda_discretization = parameter_discretization$lambda,
      init_movement_coefficients = init_movement_coefficients,
      init_stage_tx_coefficients = init_stage_tx_coefficients,
      stages_tx_from = stages_tx_from,
      stages = stages,
      population_effects = TRUE,
      validation_pct = 0.5
    )
  ),
  
  tar_target(
    name = nim_pkg_individual_val, 
    command = flatten_tags(
      tag_list = imputed_dive, 
      transition_matrices = transition_matrices, 
      template_bins = template_bins, 
      movement_types = movement_types, 
      pi_discretization = parameter_discretization$pi, 
      lambda_discretization = parameter_discretization$lambda,
      init_movement_coefficients = init_movement_coefficients,
      init_stage_tx_coefficients = init_stage_tx_coefficients,
      stages_tx_from = stages_tx_from,
      stages = stages, 
      population_effects = FALSE,
      validation_pct = 0.5
    ), 
    pattern = map(imputed_dive)
  ),
  
  tar_target(
    name = nim_fit_val,
    command = fit(
      nim_pkg = nim_pkg_val, 
      nsamples = 1e4, 
      nthin = 10, 
      max_batch_iter = 1e3, 
      max_batch_time = 5 * 60, 
      max_batch_mem = 1024^3/2, 
      sample_dir = mcmc_sample_dir,
      structural_covariates = structural_covariates
    ),
    deployment = 'worker'
  ),
  
  tar_target(
    name = nim_fit_individual_val,
    command = fit(
      nim_pkg = nim_pkg_individual_val, 
      nsamples = 1e4, 
      nthin = 10, 
      max_batch_iter = 1e3, 
      max_batch_time = 5 * 60, 
      max_batch_mem = 1024^3/2, 
      sample_dir = mcmc_sample_dir,
      structural_covariates = structural_covariates
    ),
    pattern = map(nim_pkg_individual_val), 
    deployment = 'worker'
  )
  
)
