marginalized_model_cee_prediction_script = tar_target(
  name = marginalized_model_cee_prediction,
  command = {
    
    # # load posterior samples and model package
    # tar_load(fit_marginalized_model)
    
    fit_marginalized_model = list(
      samples = "output/mcmc/fit_marginalized_model/samples.rds", 
      package = "output/mcmc/fit_marginalized_model/nim_pkg.rds"
    )
    
    samples = readRDS(fit_marginalized_model$samples)
    nim_pkg = readRDS(fit_marginalized_model$package)
    
    # set burn-in
    burn = 1:(nrow(samples)*.1)
    
    # get top-level names and groupings of variables being sampled
    sampling_targets = colnames(samples)
    sampling_target_groups = gsub(
      pattern = '\\[.*\\]',
      replacement = '',
      x = sampling_targets
    )
    sampling_groups = unique(sampling_target_groups)
    
    #
    # re-build model so that we can use it to generate intermediate objects
    #
    
    mod = nimbleModel(
      code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data,
      inits = nim_pkg$inits
    )
    
    cmod = compileNimble(mod)
    
    # verify model has a finite likelihood
    cmod$calculate()
    
    #
    # composition sample!
    #
    
    # identify posterior samples that will be used in the posterior analysis
    posterior_sample_inds = (1:nrow(samples))[-burn]
    
    # TODO: turn this into a loop over posterior samples, but maybe parallelize
    # this loop over all of the different tags since this would be a good 
    # chunk of work to parallelize over
    ind = posterior_sample_inds[1]
    
    # transfer posterior sample of model parameters to model object
    for(tgt_group in sampling_groups) {
      cmod[[tgt_group]] = samples[
        ind, sampling_targets[sampling_target_groups == tgt_group]
      ]
    }
    
    # update model components
    cmod$calculate()

    # TODO: turn this loop over subjects into an outer/primary loop
    for(subj in 1:nim_pkg$consts$n_subjects) {
      
      # get subject's id
      subject_id = nim_pkg$consts$subject_id_labels[subj]
      
      # get subject's raw data, which contains CEE information
      raw_tag = raw_sattags[[
        which(subject_id == sapply(raw_sattags, function(x) x$tag))
      ]]
      
      # get info. for the subject's final segment that was analyzed
      final_seg = nim_pkg$consts$segments[
        max(which(nim_pkg$consts$segments[,'subject_id'] == subj)),
      ]
      
      # skip processing if the animal was not clearly exposed
      if(!is.finite(raw_tag$exposure_time)) {
        err_msg = paste('There is no exposure time recorded for', 
                        subject_id,
                        'so we cannot analyze CEE response')
        warning(err_msg)
        return(list(error = err_msg, tag = subject_id))
      }
      
      # skip processing if baseline period ends considerably before exposure
      if(as.numeric(raw_tag$exposure_time) - 
         data_pkg$data$times[final_seg['end_ind']] > nim_pkg$consts$tstep) { 
        err_msg = paste('There is a large gap between the end of the baseline', 
                        'data and exposure time for', 
                        subject_id,
                        'so we cannot analyze CEE response')
        warning(err_msg)
        return(list(error = err_msg, tag = subject_id))
      }
      
      # skip processing if animal was exposed at depth
      if(data_pkg$data$covariates['depth', final_seg['end_ind']] > 
         deep_dive_depth) {
        err_msg = paste('Animal',
                        subject_id,
                        'was exposed at depth, so we skip CEE analysis')
        warning(err_msg)
        return(list(error = err_msg, tag = subject_id))
      }

      # indices of the final segment
      seg_inds = final_seg['start_ind']:final_seg['end_ind']
      
      # posterior predictive distribution for latent state before exposure
      xf = finalPred2LayerCompressed(
        obs_lik_dict = cmod$depth_tx_mat,
        obs = cmod$depth_bins[seg_inds], 
        txmat_dict = cmod$stage_tx_mat[final_seg['subject_id'], , , ],
        txmat_seq = cmod$covariateId[seg_inds],
        x0 = cmod$init_stages,
        num_obs_states = nim_pkg$consts$n_bins
      )
      
      # simulate dive until a deep depth is reached
      sim = fwd_sim_dive(
        stage = sample(
          x = 1:nim_pkg$consts$n_stages, size = 1, prob = xf
        ),
        depth_bins = nim_pkg$data$depth_bins[seg_inds], 
        covariates = data_pkg$data$covariates[, seg_inds], 
        n_timepoints = 1e3, 
        nim_pkg = nim_pkg, 
        model = cmod, 
        lon = cape_hatteras_loc['lon'], 
        lat = cape_hatteras_loc['lat'], 
        depth_threshold = deep_dive_depth,
        template_bins = template_bins,
        subj = subj
      )
      
      # extract time taken to reach a deep depth
      length(sim$depth_bins) * nim_pkg$consts$tstep
      
    }
    
    
    0
  }
)