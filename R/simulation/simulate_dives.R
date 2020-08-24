simulate_dives = function(template_bins, n_dives, obs_timesteps, 
                                 sim_pkg, pi, deep_dive_depth) {
  
  # simulate family of dives
  replicate(n_dives, {
    
    max.attempts = 1e5
    sampled.deep = FALSE
    
    attempts = 0
    while(!sampled.deep) {
      
      # sample stage transition times
      t.stages = cumsum(
        exp(rmvnorm(n = 1, mean = sim_pkg$consts$xi_prior_means[1,], 
                    sigma = sim_pkg$consts$xi_prior_covs[1,,]))
      )
      
      # simulate dive 
      d = dsdive.fwdsample.dive(
        depth.bins = template_bins, beta = pi[c(1,3)], 
        lambda = exp(sim_pkg$consts$beta_lambda_prior_mean[, 'intercept']), 
        t0 = 0, steps.max = 1e3, T1 = t.stages[1], T2 = t.stages[2]
      )
      
      # observe dive at different time intervals
      d.obs = lapply(obs_timesteps, function(tstep) {
        t.obs = seq(from = 0, to = d$times[length(d$times)] + tstep, by = tstep)
        dsdive.observe(depths = d$depths, times = d$times, stages = d$stages, 
                       t.obs = t.obs)
      })
      
      # check to see if sampled dive is deep
      attempts = attempts + 1
      sampled.deep = all(
        sapply(d.obs, function(dobs) {
          template_bins$center[max(dobs$depths)] >= deep_dive_depth
      }))
      
      # stop sampling if exceeded max samples, or got a deep dive
      if(any(attempts >= max.attempts, sampled.deep)) {
        sampled.deep = TRUE
      }
    }
    
    # label dives by their observation timestep
    names(d.obs) = obs_timesteps
    
    list(list(
      dive = d,
      dive.obs = d.obs
    ))
    
  })
   
}