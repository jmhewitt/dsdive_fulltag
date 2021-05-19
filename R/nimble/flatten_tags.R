flatten_tags = function(template_bins, lambda_discretization, stage_defs,
                        init_movement_coefficients, transition_matrices, n_pi,
                        tag_list, depth_threshold, validation_pct = 0,
                        validation_test_set = FALSE) {
  # Parameters:
  #   depth_threshold - depth used to generate prop_recent_deep covariate
  #   validation_pct - if greater than 0, then only export this percentage of 
  #     unexposed observations for model training
  #   validation_test_set - if TRUE, then export the non-validation portion of 
  #     the observations
  
  # extract dimensions
  n_bins = nrow(template_bins)
  
  # initialize flattened structures
  nim_pkg = list(
    data = list(
      constraint_data = 1,
      depths = NULL,
      covariates = NULL,
      times = NULL,
      transition_matrices = transition_matrices
    ),
    consts = list(
      lambda_discretization = lambda_discretization,
      n_stages = nrow(stage_defs),
      n_bins = nrow(template_bins),
      n_pi = n_pi,
      segments = NULL,
      stage_defs = stage_defs,
      subject_id_labels = NULL,
      n_txmat_entries = length(transition_matrices),
      n_lambda_class = length(init_movement_coefficients$lambda),
      prior_alpha = rep(1, nrow(stage_defs))
    ),
    inits = list(
      stages = NULL,
      lambda = init_movement_coefficients$lambda
    )
  )
  
  # flatten tag data
  for(tag_ind in 1:length(tag_list)) {
    
    # unwrap tag
    tag = tag_list[[tag_ind]]
    
    # store tag name
    nim_pkg$consts$subject_id_labels = c(
      nim_pkg$consts$subject_id_labels, tag$tag
    )
    
    # identify all of the observations to analyze
    tag_segments = rle(tag$gap_after)
    valid_segments = which(tag_segments$values == FALSE)
    
    # first indices of observations to analyze
    segment_starts = cumsum(c(1, tag_segments$lengths))
    
    # last index of baseline period
    baseline_end_ind = max(which(tag$times < tag$baseline_end))
    
    # last pre-exposure index
    last_pre_exposure_ind = max(which(tag$times < tag$exposure_time))
      
    # flatten data
    for(seg_ind in valid_segments) {
      # next available index in nimble package
      flat_ind = length(nim_pkg$data$depths) + 1
      # start index of data segment
      start_ind = segment_starts[seg_ind]
      # skip segment if it begins after the baseline period
      if(tag$times[start_ind] >= tag$baseline_end) {
        next
      }
      # only analyze baseline portions of tag
      end_ind = min(segment_starts[seg_ind + 1] - 1, baseline_end_ind)
      # indices to analyze
      seg_inds = seq(from = start_ind, to = end_ind, by = 1)
      # modify indices to analyze if being used for validation purposes
      if(validation_pct > 0) {
        # portion of dataset to be used for training
        train_subset = 1:ceiling(validation_pct * length(seg_inds))
        # modify data to export
        if(validation_test_set) {
          # if using as a testing dataset, only export the second portion
          seg_inds = seg_inds[-train_subset]
        } else {
          # if using as a validation dataset, only train using the first portion
          seg_inds = seg_inds[train_subset]
        }
        # update counters
        start_ind = seg_inds[1]
        end_ind = tail(seg_inds, 1)
      }
      # skip segment if there are no transitions to analyze
      if(length(seg_inds) == 1) {
        next
      }
      # copy data to package
      nim_pkg$data$depths = c(nim_pkg$data$depths, tag$depth.bin[seg_inds])
      nim_pkg$data$times = c(nim_pkg$data$times, tag$times[seg_inds])
      nim_pkg$inits$stages = c(nim_pkg$inits$stages, rep(1, length(seg_inds)))
      # proportion of recent observations at depth, as a covariate
      prop_recent_deep = c(0, sapply(2:length(seg_inds), function(i) {
        # data index to work with
        ind = seg_inds[i]
        # window at which recent observations begins
        window_start = tag$times[ind] - duration(1, units = 'hours')
        # data indices of recent observations
        past_inds = seg_inds[1:(i-1)]
        window_inds = past_inds[window_start <= tag$times[past_inds]]
        # compute proportion of recent observations spent below a depth
        sum(tag$depths[window_inds] >= depth_threshold) / length(window_inds)
      }))
      # build and add covariates
      covariates = rbind(
        intercept = rep(1, length(seg_inds)),
        daytime = tag$daytime[seg_inds],
        moonlit = tag$moonlit[seg_inds],
        prop_recent_deep = prop_recent_deep
      )
      nim_pkg$data$covariates = cbind(
        nim_pkg$data$covariates,
        covariates
      )
      # save segment information
      nim_pkg$consts$segments = rbind(
        nim_pkg$consts$segments,
        c(start_ind = flat_ind, 
          length = end_ind - start_ind + 1, 
          subject_id = tag_ind,
          end_ind = flat_ind + end_ind - start_ind)
      )
    }
  }
  
  # compute size, etc. constants
  nim_pkg$consts$n_covariates = nrow(nim_pkg$data$covariates)
  nim_pkg$consts$n_segments = nrow(nim_pkg$consts$segments)
  nim_pkg$consts$n_subjects = length(unique(
    nim_pkg$consts$segments[, 'subject_id']
  ))
  
  # discretize initial speed parameters
  for(i in 1:length(nim_pkg$inits$lambda)) {
    nim_pkg$inits$lambda_ind[i] = closestIndex(
      value = nim_pkg$inits$lambda[i],
      minval = lambda_discretization[1],
      stepsize = lambda_discretization[2],
      nvals = lambda_discretization[3]
    )
  }
  
  # initialize stage transition covariates
  nim_pkg$inits$beta_tx = array(0, dim = c(nim_pkg$consts$n_covariates,
                                           nim_pkg$consts$n_stages,
                                           nim_pkg$consts$n_stages - 1))
  
  # initialize stage transition covariate priors
  nim_pkg$consts$beta_tx_mu = numeric(nim_pkg$consts$n_stages - 1)
  nim_pkg$consts$beta_tx_cov = 1e2 * diag(nim_pkg$consts$n_stages - 1)

  nim_pkg
}
