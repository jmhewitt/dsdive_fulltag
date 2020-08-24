posterior_predictions = function(nim_pkg, sample_file, template_bins, inds, 
                                 tag, sattag_timestep) {

  0
  
  if(file.exists('nim_pkg_0.5.rds')) {
    nim_pkg = readRDS('nim_pkg_0.5.rds')
  }

  # tags with validation dives
  # tag_ids = unique(nim_pkg$consts$dive_relations_validation[, 'tag'])
  
  # extract number of train/validation dives per tag
  train_dives_per_tag = table(nim_pkg$consts$dive_relations[, 'tag'])
  test_dives_per_tag = table(nim_pkg$consts$dive_relations_validation[, 'tag'])
  
  # load posterior samples
  load(sample_file)
  
  # sample dives from posterior predictive distribution
  lapply(inds, function(ind) {
    
    #
    # simulate dive
    #
    
    # extract parameters 
    beta = samples[ind, c('pi[1]', 'pi[3]')]
    lambda = samples[ind, paste('lambda[', tag, ', ', 1:3, ']', sep = '')]
    
    # sample a traing dive from which to use stage transition times
    dive_id = sample(
      x = which(nim_pkg$consts$dive_relations[, 'tag'] == tag), 
      size = 1
    )
    
    # extract stage transition times
    stage_times = exp(
      samples[ind, paste('log_xi[', dive_id, ', ', 1:2, ']', sep ='')]
    )
    
    # sample dive
    d = dsdive.fwdsample.dive(depth.bins = template_bins, beta = beta, 
                              lambda = lambda, t0 = 0, steps.max = 1e3, 
                              T1 = stage_times[1], T2 = stage_times[2])
    
    # find dive start/end endpoints associated with dive
    endpoint_ids = nim_pkg$consts$dive_relations[
      dive_id, c('T0_endpoint', 'T3_endpoint')
    ]
    
    # extract nominal (i.e., "best guess") times for dive start/end
    obs_endpoints = apply(
      nim_pkg$consts$endpoint_priors[endpoint_ids, ], 1, mean
    )
    
    # extract sampled times for dive start/end
    est_endpoints = samples[
      ind, paste('endpoints[', endpoint_ids, ']', sep = '')
    ]
    
    # compute offsets
    offsets = obs_endpoints - est_endpoints
    names(offsets) = c('dive_start', 'dive_end')
    
    #
    # observe dive
    #
    
    # exact duration of dive
    duration = d$times[length(d$times)] - d$times[1]
    
    # build sequence of observation times
    t.obs = seq(from = max(-offsets['dive_start'], 0), 
                to = duration - offsets['dive_end'], 
                by = sattag_timestep)
    
    # observe dive
    obs = dsdive.observe(depths = d$depths, times = d$times,
                         stages = d$stages, t.obs = t.obs)
    
    # relabel observation times
    obs$times = seq(from = 0, by = sattag_timestep, 
                    length.out = length(t.obs))
    
    # compute observed stage durations
    stages.dur = diff(c(0, obs$times[c(FALSE, diff(obs$stages)==1)],
                        obs$times[length(obs$times)]))
    
    # package results
    list(
      dive = d,
      dive.obs = obs,
      stages.dur = stages.dur,
      offsets = offsets
    )
  })

}