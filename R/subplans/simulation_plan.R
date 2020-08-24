simulation_plan = drake_plan(
  
  # simulate dives
  sim_dives = target(
    simulate_dives(
      template_bins = template_bins, 
      n_dives = 100, 
      obs_timesteps = c(30, 60, 300),
      sim_pkg = nim_pkg_0,
      pi = c(.97, .5, .03),
      deep_dive_depth = deep_dive_depth
    ), 
    seed = 2020
  ),
  
  # flatten simulated dives 
  nim_sim_data = target(
    flatten_simulation_data(sim_dives, template_bins, timestep),
    transform = map(timestep = !!(c(30,60,300)))
  ),
  
  # build priors and initialize parameters
  nim_sim_pkg = target(
    priors_and_inits(nim_sim_data, template_bins, timestep, expm_delta),
    transform = map(nim_sim_data, 
                    timestep = !!(c(30,60,300)), 
                    .names = !!(paste('nim_sim_pkg', c(30,60,300), sep = '_'))),
    trigger = trigger(condition = FALSE, mode = 'blacklist')
  ),
  
  # fit models
  mcmc_samples_sim = target(
    fit(nim_sim_pkg, mcmc_sample_dir, niter, ncheckpoints, .id_chr, 
        default_ranef_samplers = TRUE),
    transform = map(nim_sim_pkg, 
                    .names = !!(paste('mcmc_samples_sim', c(30,60,300), 
                                      sep = '_'))),
    format = 'file',
    trigger = trigger(condition = FALSE, mode = 'blacklist')
  ),
  
  # summarize simulation results
  sim_summary = target(
    sim_summary(sim_pkg = nim_pkg_0, pi = c(.97, .5, .03), 
                sample_file = mcmc_samples_sim,
                nburn = nburn),
    transform = map(mcmc_samples_sim, 
                    .names = !!(paste('sim_summary', c(30,60,300), 
                                      sep = '_')))
  )
  
)
