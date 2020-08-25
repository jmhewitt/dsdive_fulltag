nimble_plan = drake_plan(
  
  # location to store posterior samples
  mcmc_sample_dir = file_out(!!file.path('output', 'mcmc')),
  
  # constant to add to infinitesimal generator diagonals for numerical speedup
  expm_delta = 1e-10,
  
  # number of posterior samples to generate
  niter = 1e4,
  
  # number of posterior samples to discard
  nburn = 5e3,
  
  # number of times to save partial output while sampling
  ncheckpoints = 20,
  
  # reformat dive data for analysis in nimble (full dataset and validaiton set)
  nim_data = target(
    flatten_tag_data(depth_files[1], dive_endpoints[1], template_bins, 
                     tag_sex, cee_starts, validation_proportion = 0, 
                     mcmc_sample_dir), 
    format = 'file'
  ),
  
  # build priors and initialize parameters
  nim_pkg = target(
    priors_and_inits(nim_data, template_bins, sattag_timestep, expm_delta, 
                     mcmc_sample_dir),
    format = 'file'
  ),

  # fit models
  mcmc_samples = target(
    fit(nim_pkg, mcmc_sample_dir, niter, ncheckpoints),
    format = 'file'
  )
  
)

