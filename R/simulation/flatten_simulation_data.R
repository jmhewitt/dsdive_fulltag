flatten_simulation_data = function(sim_dives, template_bins, timestep) {
  
  # get number of simulation dives
  n_dives = length(sim_dives)
  
  # initialize flattened structures
  nim_pkg = list(
    data = list(
      depths = NULL,
      times = NULL
    ),
    consts = list(
      endpoint_priors = matrix(nrow = n_dives * 2, ncol = 2),
      dive_relations = matrix(nrow = n_dives, ncol = 5),
      dive_relations_validation = NULL,
      tag_covariates = NULL
    )
  )
  
  colnames(nim_pkg$consts$endpoint_priors) = c('t_lwr', 't_upr')
  colnames(nim_pkg$consts$dive_relations) = c(
    'T0_endpoint', 'T3_endpoint', 'tag', 'depth_first', 'depth_last'
  )
  
  # "intercept" sex covariate for the record
  nim_pkg$consts$tag_covariates = c(0, 0)
  
  # flatten each simulated dive
  for(dive.id in 1:length(sim_dives)) {
    
    # extract dive observations
    d = sim_dives[[dive.id]]$dive.obs[[as.character(timestep)]]
    
    # construct ids for endpoints
    flat_endpoint_ids = 2 * dive.id + c(-1, 0)
    
    # set priors for dive start endpoint
    nim_pkg$consts$endpoint_priors[flat_endpoint_ids[1], ] = 
      timestep * c(-1,1) + d$times[1]
    
    # set priors for dive end endpoint
    nim_pkg$consts$endpoint_priors[flat_endpoint_ids[2], ] = 
      timestep * c(-1,1) + d$times[length(d$times)]
    
    # get index at which first depth will be stored
    flat_depth_ind = length(nim_pkg$data$depths) + 1
    
    # merge observed depths and times
    nim_pkg$data$depths = c(nim_pkg$data$depths, d$depths)
    nim_pkg$data$times = c(nim_pkg$data$times, d$times)
    
    # merge dive information in flattened structure
    nim_pkg$consts$dive_relations[
      dive.id, c('T0_endpoint', 'T3_endpoint', 'tag', 'depth_first', 
                 'depth_last')
    ] = c(flat_endpoint_ids, 1, flat_depth_ind, length(nim_pkg$data$depths))
    
  }
  
  # extract sizes
  nim_pkg$consts$N_tags = 1
  nim_pkg$consts$N_dives = nrow(nim_pkg$consts$dive_relations)
  nim_pkg$consts$N_endpoints = nrow(nim_pkg$consts$endpoint_priors)
  nim_pkg$consts$N_bins = nrow(template_bins)
  
  # export bin widths
  nim_pkg$consts$widths = template_bins$halfwidth * 2
  
  # remove validation data placeholder
  nim_pkg$consts$dive_relations_validation = NULL
  
  nim_pkg
}
