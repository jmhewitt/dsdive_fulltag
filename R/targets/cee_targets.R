cee_targets = list(
  
  tar_target(
    cee_dive_response_targets, 
    c("ZcTag093", "ZcTag095", "ZcTag096", "ZcTag097")
  ),
  
  tar_target(
    name = cee_dive_response_probs, 
    command = {
      
      # location of MCMC files (also used as output path)
      path = nim_fit
      
      # posterior parameter sample files
      param_sample_files = dir(
        path = path, pattern = 'parameter_samples_[0-9]+', full.names = TRUE
      )
      
      # posterior stage sample files
      stage_sample_files = dir(
        path = path, pattern = 'stage_samples_[0-9]+', full.names = TRUE
      )
      
      # column labels for posterior parameter samples
      param_label_file = dir(
        path = path, pattern = 'parameter_samples_column', full.names = TRUE
      )
      
      # load posterior parameter samples 
      param_samples = do.call(rbind, lapply(param_sample_files, function(f) {
        readRDS(f)
      }))
      
      # load posterior stage samples
      stage_samples = do.call(rbind, lapply(stage_sample_files, function(f) {
        readRDS(f)
      }))
      
      # load data package
      nim_pkg = readRDS(
        file = dir(path = path, pattern = 'nim_pkg', full.names = TRUE)
      )
      
      # label posterior parameter samples
      colnames(param_samples) = readRDS(param_label_file)
      
      burn = 1:1e3
      
      # indices of posterior samples to process
      post_samples = (1:nrow(param_samples))[-burn]
      
      # match tag id in model
      subj_id = which(
        cee_dive_response_targets == nim_pkg$consts$subject_id_labels
      )
      
      # extract pre-exposure data index
      pre_exposure_ind = data.frame(nim_pkg$consts$segments) %>% 
        filter(subject_id == subj_id) %>% 
        select(end_ind) %>% 
        unlist() %>% 
        max()
      
      # subset posterior samples according to parallelization task
      post_samples = post_samples
      
      # draw posterior predictive samples of time to next deep observation
      post_pred_samples = sapply(post_samples, function(sample_ind) {
          
          # extract stage transition coefficients
          betas_tx = array(dim = c(nim_pkg$consts$n_covariates,
                                   nim_pkg$consts$n_stages,
                                   nim_pkg$consts$n_stages - 1))
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
          
          # forward-simulate dive from predictive distribution
          init_inds = seq(to = pre_exposure_ind, by = 1, length.out = 12)
          fwd_sim = fwd_sim_to_deep(
            stages = stage_samples[sample_ind, init_inds], 
            depths = nim_pkg$data$depths[init_inds], 
            covariates = nim_pkg$data$covariates[, init_inds], 
            n_max = 1e3, 
            nim_pkg = nim_pkg,
            lambda = param_samples[sample_ind, c('lambda[1]', 'lambda[2]')], 
            betas_tx = betas_tx, 
            template_bins = template_bins, 
            times = as.POSIXct(
              x = nim_pkg$data$times[init_inds],
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
      
      list(list(
        tag = cee_dive_response_targets,
        samples = post_pred_samples
      ))
      
    }, 
    pattern = map(cee_dive_response_targets)
  )
  
  
)
