nimble_targets = list(
  
  # merge data for analysis into a single structure that nimble can understand
  tar_target(
    name = data_pkg, 
    command = flatten_tags(
      template_bins = template_bins,
      tag_list = raw_sattags,
      depth_threshold = deep_dive_depth,
      sattag_timestep = sattag_timestep,
      repeated_surface_bin_break = 3,
      min_segment_length = 60
    )
  ),
  
  # define discrete movement classes/states, in which each class is associated
  # with one of the modeled lambda values, and one of the discretized pi values
  tar_target(
    name = movement_classes, 
    command = {
      res = list(
        # define how we should interpret the parameters
        lambda_class = c('slow', 'fast'),
        pi_class = c('descent', 'free', 'ascent')
      )
      
      # define how the parameters should be combined to form a movement class
      res$stage_defs = rbind(
        c(pi_ind = which(res$pi_class == 'descent'), 
          lambda_ind = which(res$lambda_class == 'slow')),
        c(pi_ind = which(res$pi_class == 'descent'), 
          lambda_ind = which(res$lambda_class == 'fast')),
        c(pi_ind = which(res$pi_class == 'free'), 
          lambda_ind = which(res$lambda_class == 'slow')),
        c(pi_ind = which(res$pi_class == 'ascent'), 
          lambda_ind = which(res$lambda_class == 'slow')),
        c(pi_ind = which(res$pi_class == 'ascent'), 
          lambda_ind = which(res$lambda_class == 'fast'))
      )
      
      # generate labels for the movement classes
      rownames(res$stage_defs) = apply(res$stage_defs, 1, function(r) {
        paste(res$lambda_class[r['lambda_ind']],
              res$pi_class[r['pi_ind']],
              sep = '_'
        )
      })
      
      res
    }
  ),
  
  tar_target(
    name = covariate_tx_control,
    command = list(
      deep_depth = deep_dive_depth,
      window_len = 3600,
      obs_freq = sattag_timestep,
      poly_degree = 3,
      vertical_scale = 500,
      vertical_center = 1650
    )
  ),
  
  fit_marginalized_model_script,
  
  marginalized_model_posterior_diagnostic_script,
  
  marginalized_model_cee_prediction_script,
  
  parameter_interpretation_pattern_script
  
)
