flatten_tags = function(tag_list, transition_matrices, movement_types,
                        pi_discretization, lambda_discretization, 
                        template_bins, init_movement_coefficients,
                        init_stage_tx_coefficients, stages_tx_from, stages,
                        population_effects, validation_pct = 0) {
  
  # extract dimensions
  n_bins = nrow(template_bins)
  
  # extract initial coefficients
  alpha = init_movement_coefficients$alpha
  beta = init_movement_coefficients$beta
  betas_tx = init_stage_tx_coefficients
  
  # initialize flattened structures
  nim_pkg = list(
    data = list(
      depths = NULL,
      times = NULL,
      stages = NULL,
      stage_supports = NULL,
      covariates = NULL,
      surface_bin = NULL,
      transition_matrices = transition_matrices
    ),
    consts = list(
      segments = NULL,
      n_bins = n_bins,
      movement_types = movement_types,
      pi_discretization = pi_discretization,
      lambda_discretization = lambda_discretization,
      n_pi = as.integer(pi_discretization[, 'nvals']),
      n_lambda = as.integer(lambda_discretization[, 'nvals']),
      subject_id_labels = NULL,
      betas_tx_stage_from = as.integer(stages_tx_from),
      population_effects = population_effects
    ),
    inits = list(
      # population-level effects
      alpha_mu = alpha,
      beta_mu = beta,
      betas_tx_mu = betas_tx,
      alpha_var = matrix(1, nrow = nrow(alpha), ncol = ncol(alpha)),
      beta_var = matrix(1, nrow = nrow(beta), ncol = ncol(beta)),
      betas_tx_var = matrix(1, nrow = nrow(betas_tx), ncol = ncol(betas_tx))
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
    tag_segments = rle(tag$data_gaps)
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
      # if using as a validation dataset, only train using the first portion
      if(validation_pct > 0) {
        train_subset = 1:ceiling(validation_pct * length(seg_inds))
        seg_inds = seg_inds[train_subset]
        end_ind = tail(seg_inds, 1)
      }
      # skip segment if there are no transitions to analyze
      if(length(seg_inds) == 1) {
        next
      }
      # copy data to package
      nim_pkg$data$depths = c(nim_pkg$data$depths, tag$depth.bin[seg_inds])
      nim_pkg$data$stages = c(nim_pkg$data$stages, tag$stages[seg_inds])
      nim_pkg$data$stage_supports = cbind(
        nim_pkg$data$stage_supports,
        tag$stage_support[,seg_inds]
      )
      nim_pkg$data$times = c(nim_pkg$data$times, tag$times[seg_inds])
      # build and add covariates
      covariates = rbind(
        intercept = rep(1, length(seg_inds)),
        depth =  tag$depths[seg_inds],
        deep_depth = tag$depths[seg_inds] >= 800,
        shallow_depth = tag$depths[seg_inds] < 800,
        non_surface_bin = tag$depth.bin[seg_inds] > 1,
        surface_bin = tag$depth.bin[seg_inds] == 1,
        time_since_surface = rep(0, length(seg_inds)),
        all_shallow_depths_since_surface = rep(0, length(seg_inds)),
        daytime = tag$daytime[seg_inds],
        moonlit = tag$moonlit[seg_inds]
      )
      covariates['time_since_surface',] =  (
        1 + time_in_state(covariates['non_surface_bin',], ncol(covariates))
      ) * covariates['non_surface_bin',]
      all_shallow = TRUE
      for(i in 1:ncol(covariates)) {
        if(covariates['surface_bin',i] == TRUE) {
          all_shallow = TRUE
        } else if(covariates['deep_depth',i] == TRUE) {
          all_shallow = FALSE
        }
        covariates['all_shallow_depths_since_surface',i] = all_shallow
      }
      nim_pkg$data$covariates = cbind(
        nim_pkg$data$covariates,
        covariates
      )
      nim_pkg$data$surface_bin = c(
        nim_pkg$data$surface_bin,
        tag$depth.bin[seg_inds] == 1
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
  nim_pkg$consts$n_timepoints = length(nim_pkg$data$depths)
  nim_pkg$consts$n_stage_txs = ncol(nim_pkg$inits$betas_tx_mu)
  nim_pkg$consts$n_stages = length(stages)
  nim_pkg$consts$n_txmat_entries = length(nim_pkg$data$transition_matrices)
  nim_pkg$consts$n_txmat_types = nrow(nim_pkg$consts$pi_discretization)
  nim_pkg$consts$n_covariates = nrow(nim_pkg$data$covariates)
  nim_pkg$consts$n_segments = nrow(nim_pkg$consts$segments)
  nim_pkg$consts$n_subjects = length(unique(
    nim_pkg$consts$segments[, 'subject_id']
  ))
  
  # add initial individual-level random effects
  nim_pkg$inits$alpha = array(
    data = alpha, 
    dim = c(nrow(alpha), ncol(alpha), nim_pkg$consts$n_subjects)
  )
  rownames(nim_pkg$inits$alpha) = rownames(alpha)
  colnames(nim_pkg$inits$alpha) = colnames(alpha)
  nim_pkg$inits$beta = array(
    data = beta, 
    dim = c(nrow(beta), ncol(beta), nim_pkg$consts$n_subjects)
  )
  rownames(nim_pkg$inits$beta) = rownames(beta)
  colnames(nim_pkg$inits$beta) = colnames(beta)
  nim_pkg$inits$betas_tx = array(
    data = betas_tx, 
    dim = c(nrow(betas_tx), ncol(betas_tx), nim_pkg$consts$n_subjects)
  )
  rownames(nim_pkg$inits$betas_tx) = rownames(betas_tx)
  colnames(nim_pkg$inits$betas_tx) = colnames(betas_tx)
  
  # # initial group-level intercept contribution for stage transition params.
  # nim_pkg$inits$betas_tx_stage_offset = t(apply(betas_tx, 1, function(x) {
  #   aggregate(x, list(nim_pkg$consts$betas_tx_stage_from), mean)[,'x']
  # }))
  # colnames(nim_pkg$inits$betas_tx_stage_offset) = names(stages)
  # rownames(nim_pkg$inits$betas_tx_stage_offset) = rownames(
  #   nim_pkg$data$covariates
  # )
  
  # # adjust means for stage offsets/intercepts
  # nim_pkg$inits$betas_tx_mu = 
  #   nim_pkg$inits$betas_tx_stage_offset[,nim_pkg$consts$betas_tx_stage_from] - 
  #   nim_pkg$inits$betas_tx_mu
  # colnames(nim_pkg$inits$betas_tx_mu) = colnames(betas_tx)
  # 
  # nim_pkg$inits$betas_tx_eps = array(
  #   data = nim_pkg$inits$betas_tx_mu, 
  #   dim = c(nrow(nim_pkg$inits$betas_tx_mu), ncol(nim_pkg$inits$betas_tx_mu), 
  #           nim_pkg$consts$n_subjects)
  # )
  # rownames(nim_pkg$inits$betas_tx_eps) = rownames(nim_pkg$inits$betas_tx_mu)
  # colnames(nim_pkg$inits$betas_tx_eps) = colnames(nim_pkg$inits$betas_tx_mu)
  
  # add information about key stages and covariates
  nim_pkg$consts$intercept_covariate = which(
    rownames(nim_pkg$data$covariates) == 'intercept'
  )
  nim_pkg$consts$ascent_like_stages = stages[
    c('deep_ascent', 'shallow_ascent', 'free_surface')
  ]
  nim_pkg$consts$n_ascent_like_stages = length(
    nim_pkg$consts$ascent_like_stages
  )
  
  # initial parameters become priors if population-level effects aren't est'd.
  if(population_effects == FALSE) {
    nim_pkg$inits$betas_tx_mu = 0 * nim_pkg$inits$betas_tx_mu
    nim_pkg$inits$beta_mu = 0 * nim_pkg$inits$beta_mu
    nim_pkg$inits$beta_var = matrix(1e2, nrow = nrow(beta), ncol = ncol(beta))
    nim_pkg$inits$betas_tx_var = matrix(
      1e2, nrow = nrow(betas_tx), ncol = ncol(betas_tx)
    )
  }

  nim_pkg
}
