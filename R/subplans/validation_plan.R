validation_plan = drake_plan(
  
  # location to store validation output files
  validation_output_dir = file_out(!!file.path('output', 'validation')),
  
  # posterior predictive samples to be used for validation
  posterior_sample_inds = nburn:niter,
  
  # simulate dives from posterior predictive distribution
  posterior_predictive_validation_samples = target(
    posterior_predictions(nim_pkg = nim_pkg_0.5,
                          sample_file = mcmc_samples_nim_pkg_0.5,
                          template_bins = template_bins,
                          inds = posterior_sample_inds,
                          tag = tag,
                          sattag_timestep = sattag_timestep),
    transform = map(tag = !!(1:7), .tag_out = tag)
  ),
  
  # compare validation dives to posterior predictive distribution
  validation_statistics = target(
    validation_statistics(nim_pkg = nim_pkg_0.5,
                          samples = posterior_predictive_validation_samples,
                          tag = id_chr(),
                          deep_dive_depth = deep_dive_depth,
                          template_bins = template_bins,
                          sattag_timestep = sattag_timestep,
                          validation_output_dir = validation_output_dir,
                          tag_sex = tag_sex),
    transform = map(tag, 
                    .names = !!(paste('validation_statistics', 1:7, sep = '_')))
  )
  
)

