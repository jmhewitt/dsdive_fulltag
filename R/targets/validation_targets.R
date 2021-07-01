validation_targets = list(
  
  tar_target(
    name = nim_pkg_val, 
    command = flatten_tags(
      template_bins = template_bins, 
      lambda_discretization = parameter_discretization$lambda, 
      stage_defs = movement_classes$stage_defs, 
      init_movement_coefficients = list(lambda = rep(1,3)), 
      transition_matrices = transition_matrices, 
      n_pi = length(parameter_discretization$pi), 
      tag_list = raw_sattags,
      depth_threshold = deep_dive_depth,
      validation_pct = 0.5, 
      repeated_surface_bin_break = repeated_surface_bin_break, 
      min_segment_length = min_segment_length
    )
  ),
  
  tar_target(
    name = nim_pkg_val_test, 
    command = flatten_tags(
      template_bins = template_bins, 
      lambda_discretization = parameter_discretization$lambda, 
      stage_defs = movement_classes$stage_defs, 
      init_movement_coefficients = list(lambda = rep(1,3)), 
      transition_matrices = transition_matrices, 
      n_pi = length(parameter_discretization$pi), 
      tag_list = raw_sattags,
      depth_threshold = deep_dive_depth,
      validation_pct = 0.5,
      validation_test_set = TRUE, 
      repeated_surface_bin_break = repeated_surface_bin_break, 
      min_segment_length = min_segment_length
    )
  ),
  
  tar_target(
    name = nim_fit_val,
    command = fit(
      nim_pkg = nim_pkg_val, 
      nsamples = 1e4, 
      nthin = 10, 
      max_batch_iter = 1e3, 
      max_batch_time = 30 * 60, 
      max_batch_mem = 1024^3/2, 
      sample_dir = mcmc_sample_dir
    ),
    deployment = 'worker'
  ),
  
  tar_target(
    name = validation_dives,
    command = extract_validation_dives(
      tag_list = imputed_dive,
      validation_pct = 0.5
  )),
  
  
  # batch settings for validation dives
  tar_target(validation_dive_batch_size, 10),
  tar_target(
    name = validation_batch_starts, 
    command = seq(
      from = 1, 
      to = length(validation_dives$data$dives), 
      length.out = 100
    )
  ),
  
  # tar_target(val_timepoints, seq(from = 5, to = 30, by = 5)),
  tar_target(val_timepoints, seq(from = 1, to = 15, by = 1)),
  
  tar_target(
    name = validate_dive_preds,
    command = validate_dive_predictions(
      post_output_dir = file.path('output', 'mcmc', 'nim_fit_val'),
      validation_dives = validation_dives$data$dives[
        validation_batch_starts + 0:(validation_dive_batch_size -1)
      ],
      burn = 1e3,
      n_timepoints = val_timepoints
    ),
    pattern = cross(val_timepoints, validation_batch_starts),
    deployment = 'worker'
  ),
    
  tar_target(val_dive_length_timepoints, seq(from = 5, to = 60, by = 5)),
  
  tar_target(
    name = validate_dive_length_predictions,
    command = validate_dive_length_predictions(
      post_output_dir = file.path('output', 'mcmc', 'nim_fit_val'),
      validation_dives = validation_dives$data$dives[
        validation_batch_starts + 0:(validation_dive_batch_size -1)
      ],
      burn = 1e3,
      n_timepoints = val_dive_length_timepoints,
      timestep = sattag_timestep/imputation_factor, 
      lon = cape_hatteras_loc['lon'],
      lat = cape_hatteras_loc['lat'],
      template_bins = template_bins
    ),
    pattern = cross(val_dive_length_timepoints, validation_batch_starts),
    deployment = 'worker'
  ),
  
  tar_target(
    name = validation_df_deep_surv,
    command = {

      # data indices for validation dives
      val_inds = which(
        # dives must start on surface
        (c(NA, 
           nim_pkg_val_test$data$depths[
             1:(length(nim_pkg_val_test$data$depths)-1)
           ]) == 1) &
        # ...and be at depth bin 5 at next observation
        (nim_pkg_val_test$data$depths == 5) &
        # previous-depth covariate must have "full records"
        (nim_pkg_val_test$data$n_previous_used == 12) &
        # dives must be strictly daytime
        (nim_pkg_val_test$data$covariates['daytime',] == 1) &
        (nim_pkg_val_test$data$covariates['moonlit',] == 0)
      )

      # determine time to next deepest depth, or whether response is censored
      val_responses = do.call(rbind, lapply(val_inds, function(dive_ind) {
        # determine segment in which dive was recorded
        seg_ind = which(
          (nim_pkg_val_test$consts$segments[,'start_ind'] <= dive_ind) &
            (dive_ind <= nim_pkg_val_test$consts$segments[,'end_ind'])
        )
        # data indices in which to search for deep depth observation
        test_inds = seq(
          from = dive_ind + 1,
          to = nim_pkg_val_test$consts$segments[seg_ind, 'end_ind']
        )
        # package results
        data.frame(
          # return time of first deep observation (Inf if censored)
          deep_time = min(which(
            template_bins$center[nim_pkg_val_test$data$depths[test_inds]] >= 
              deep_dive_depth  
          )),
          # num. observations checked (for partial analyses of censored data)
          nobs = length(test_inds)
        )
      }))

      # package results
      data.frame(
        pkg_data_ind = val_inds,
        val_responses,
        t(nim_pkg_val_test$data$covariates[, val_inds])
      )
    }
  ),
  
  # define parallelization for validation
  tar_target(n_validation_tasks, 100),
  tar_target(validation_tasks, 1:n_validation_tasks),
  
  tar_target(
    name = validate_deep_surv_samples,
    command = {
      
      # location of MCMC files (also used as output path)
      # path = nim_fit_val
      path = file.path('output', 'mcmc', 'nim_fit_val')
      
      # posterior parameter sample files
      param_sample_files = dir(
        path = path, pattern = 'parameter_samples_[0-9]+', full.names = TRUE
      )
      
      # column labels for posterior parameter samples
      param_label_file = dir(
        path = path, pattern = 'parameter_samples_column', full.names = TRUE
      )
      
      # load posterior parameter samples 
      param_samples = do.call(rbind, lapply(param_sample_files, function(f) {
        readRDS(f)
      }))
      
      # label posterior parameter samples
      colnames(param_samples) = readRDS(param_label_file)
      
      burn = 1:1e3
      
      # indices of posterior samples to process
      post_samples = (1:nrow(param_samples))[-burn]
      
      # select posterior samples to process based on parallelization task num.
      
      # determine parallelization
      parallel_tasks = parallel::splitIndices(
        nx = length(post_samples), 
        ncl = n_validation_tasks
      )
      
      # subset posterior samples according to parallelization task
      post_samples = post_samples[parallel_tasks[[validation_tasks]]]

      # "load" validation dataset
      df.val = validation_df_deep_surv
      
      #
      # draw posterior validation samples
      #
      
      validation_samples = do.call(
        rbind, lapply(post_samples, function(sample_ind) {
        
          # extract stage transition coefficients
          betas_tx = array(dim = c(nim_pkg_val_test$consts$n_covariates,
                                   nim_pkg_val_test$consts$n_stages,
                                   nim_pkg_val_test$consts$n_stages - 1))
          for(i in 1:nrow(betas_tx)) {
            for(j in 1:ncol(betas_tx)) {
              for(k in 1:dim(betas_tx)[3]) {
                betas_tx[i,j,k] = param_samples[
                  sample_ind, paste('beta_tx[', i, ', ', j, ', ', k, ']', 
                                    sep = '')
                ]
              }
            }
          }
          
          # sample posterior predictive distribution for each validation dive
          pred_samples = sapply(df.val$pkg_data_ind, function(dive_ind) {
            
            # selected validation dives are assumed to be slow descents b/c
            # of m/s speed implied by depth bin 1 -> 5 transition
            stage_start = which(
              rownames(movement_classes$stage_defs) == 'slow_descent'
            )
            
            # forward-simulate dive from predictive distribution
            init_inds = seq(to = dive_ind, by = 1, length.out = 12)
            fwd_sim = fwd_sim_to_depth(
              stages = rep(stage_start, 12), 
              depths = nim_pkg_val_test$data$depths[init_inds], 
              covariates = nim_pkg_val_test$data$covariates[, init_inds], 
              n_max = 1e3, 
              nim_pkg = nim_pkg_val_test,
              lambda = param_samples[sample_ind, c('lambda[1]', 'lambda[2]')], 
              betas_tx = betas_tx, 
              template_bins = template_bins, 
              times = as.POSIXct(
                x = nim_pkg_val_test$data$times[init_inds],
                tz = 'UTC', 
                origin = '1970-01-01 00:00.00 UTC'
              ), 
              timestep = sattag_timestep, 
              lon = cape_hatteras_loc['lon'], lat = cape_hatteras_loc['lat'], 
              depth_threshold = deep_dive_depth
            )
          
            # sampled time to deep depth  
            length(fwd_sim$depths) - 1
          })
          
          # return samples
          pred_samples
      }))
      
      # package results
      res = list(
        samples = validation_samples,
        posterior_sample_inds = post_samples
      )
      
      # save results
      saveRDS(res, file = file.path(path, paste(tar_name(), '.rds', sep = '')))
      
      # return output directory
      path
    }, 
    pattern = map(validation_tasks), 
    deployment = 'worker'
  ),
  
  tar_target(
    name = validate_deep_surv_distn, 
    command = {
      
      # directory where posterior samples exist
      path = validate_deep_surv_samples
      # path ='output/mcmc/nim_fit_val/'
      
      # list of validation files
      flist = dir(path = path, pattern = 'validate', full.names = TRUE)
      
      # summarize posterior predictive samples with frequency tables
      aggregated_subsets = lapply(flist, function(f) {
        apply(readRDS(f)$samples, 2, function(x) {
          res = table(x)
          data.frame(x = as.numeric(names(res)), Freq = as.numeric(res))
        })
      })
      
      # convenience function to assist in aggregating frequency tables
      merge_tables = function(x1, x2) {
        # merge tables to get complete list of index values
        m = x1 %>% full_join(x2, by = 'x')
        # set na's in frequency counts to 0
        m[is.na(m)] = 0
        # return aggregated counts
        m %>%
          mutate(Freq = Freq.x + Freq.y) %>% 
          select(x, Freq)
      }
      
      # extract posterior predictive distributions
      post_cdfs = lapply(1:length(aggregated_subsets[[1]]), function(ind) {
        # merge subsets of posterior samples
        pmf = aggregated_subsets[[1]][[ind]]
        if(length(aggregated_subsets) > 1) {
          for(i in 2:length(aggregated_subsets)) {
            pmf = merge_tables(x1 = pmf, x2 = aggregated_subsets[[i]][[ind]])
          }
        }
        # normalize and sort pmf
        pmf$Freq = pmf$Freq / sum(pmf$Freq)
        pmf = pmf[order(pmf$x),]
        # return result
        pmf
      })
      
      # return result
      post_cdfs
    })
  
)
