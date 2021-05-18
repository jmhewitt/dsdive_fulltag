validation_targets = list(
  
  tar_target(
    name = nim_pkg_val, 
    command = flatten_tags(
      template_bins = template_bins, 
      lambda_discretization = parameter_discretization$lambda, 
      stage_defs = movement_classes$stage_defs, 
      init_movement_coefficients = list(lambda = rep(1,3)), 
      transition_matrices = transition_matrices, 
      n_pi = length(parameter_discretization$pi), 
      tag_list = raw_sattags,
      depth_threshold = deep_dive_depth,
      validation_pct = 0.5
    )
  ),
  
  tar_target(
    name = nim_pkg_val_test, 
    command = flatten_tags(
      template_bins = template_bins, 
      lambda_discretization = parameter_discretization$lambda, 
      stage_defs = movement_classes$stage_defs, 
      init_movement_coefficients = list(lambda = rep(1,3)), 
      transition_matrices = transition_matrices, 
      n_pi = length(parameter_discretization$pi), 
      tag_list = raw_sattags,
      depth_threshold = deep_dive_depth,
      validation_pct = 0.5,
      validation_test_set = TRUE
    )
  ),
  
  tar_target(
    name = nim_fit_val,
    command = fit(
      nim_pkg = nim_pkg_val, 
      nsamples = 1e4, 
      nthin = 10, 
      max_batch_iter = 1e3, 
      max_batch_time = 30 * 60, 
      max_batch_mem = 1024^3/2, 
      sample_dir = mcmc_sample_dir
    ),
    deployment = 'worker'
  ),
  
  tar_target(
    name = validation_dives,
    command = extract_validation_dives(
      tag_list = imputed_dive,
      validation_pct = 0.5
  )),
  
  
  # batch settings for validation dives
  tar_target(validation_dive_batch_size, 10),
  tar_target(
    name = validation_batch_starts, 
    command = seq(
      from = 1, 
      to = length(validation_dives$data$dives), 
      length.out = 100
    )
  ),
  
  # tar_target(val_timepoints, seq(from = 5, to = 30, by = 5)),
  tar_target(val_timepoints, seq(from = 1, to = 15, by = 1)),
  
  tar_target(
    name = validate_dive_preds,
    command = validate_dive_predictions(
      post_output_dir = file.path('output', 'mcmc', 'nim_fit_val'),
      validation_dives = validation_dives$data$dives[
        validation_batch_starts + 0:(validation_dive_batch_size -1)
      ],
      burn = 1e3,
      n_timepoints = val_timepoints
    ),
    pattern = cross(val_timepoints, validation_batch_starts),
    deployment = 'worker'
  ),
    
  tar_target(val_dive_length_timepoints, seq(from = 5, to = 60, by = 5)),
  
  tar_target(
    name = validate_dive_length_predictions,
    command = validate_dive_length_predictions(
      post_output_dir = file.path('output', 'mcmc', 'nim_fit_val'),
      validation_dives = validation_dives$data$dives[
        validation_batch_starts + 0:(validation_dive_batch_size -1)
      ],
      burn = 1e3,
      n_timepoints = val_dive_length_timepoints,
      timestep = sattag_timestep/imputation_factor, 
      lon = cape_hatteras_loc['lon'],
      lat = cape_hatteras_loc['lat'],
      template_bins = template_bins
    ),
    pattern = cross(val_dive_length_timepoints, validation_batch_starts),
    deployment = 'worker'
  )
  
)
