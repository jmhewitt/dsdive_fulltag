cee_targets = list(
  
  tar_target(
    cee_tag_targets, 
    c("ZcTag069", "ZcTag085", "ZcTag086", "ZcTag093", 
      "ZcTag095", "ZcTag096", "ZcTag097")
  ),
  
  tar_target(
    name = dive_length_predictions,
    command = predict_dive_length(
      post_output_dir = file.path('output', 'mcmc', 'nim_fit'),
      tag = cee_tag_targets,
      burn = 1:1e3,
      timestep = sattag_timestep/imputation_factor,
      lon = cape_hatteras_loc['lon'],
      lat = cape_hatteras_loc['lat'],
      template_bins = template_bins
    ),
    pattern = map(cee_tag_targets),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  )
  
)
