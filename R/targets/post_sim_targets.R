post_sim_targets = list(

  tar_target(
    name = covariate_combinations,
    command = expand.grid(
      prop_recent_deep = (1:12)/12,
      daytime = c(0,1),
      moonlit = c(0,1)
    )
  ),
  
  tar_target(
    name = deep_dive_time_preds_wrt_covariates, 
    command = {
      
      # location of MCMC files (also used as output path)
      # path = nim_fit
      path = 'output/mcmc/nim_fit/'
      
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
      
      # load data package
      nim_pkg = readRDS(
        file = dir(path = path, pattern = 'nim_pkg', full.names = TRUE)
      )
      
      # label posterior parameter samples
      colnames(param_samples) = readRDS(param_label_file)
      
      burn = 1:2e3
      
      # indices of posterior samples to process
      post_samples = (1:nrow(param_samples))[-burn]
      
      # subset posterior samples according to parallelization task
      post_samples = post_samples
      
      # build initial covariate matrix
      covariates = rbind(
        intercept = rep(1, 12),
        daytime = rep(covariate_combinations$daytime, 12), 
        moonlit = rep(covariate_combinations$moonlit, 12),
        prop_recent_deep = pmax(
          covariate_combinations$prop_recent_deep - (0:11)/12, 
          0
        ),
        prop_recent_deep3 = rep(0,12)
      )
      covariates['prop_recent_deep3',] = (
        covariates['prop_recent_deep',] - 0.5
      )^3
      
      # build initial depth vector
      depths = rep(x = c(1,16), 
          c(12 - covariate_combinations$prop_recent_deep * 12,
            covariate_combinations$prop_recent_deep * 12)
      )
      
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
        fwd_sim = fwd_sim_to_depth_fixed_covs(
          stages = rep(
            which(rownames(movement_classes$stage_defs) == 'slow_descent'), 
            12
          ),
          depths = depths, 
          covariates = covariates, 
          n_max = 1e3, 
          nim_pkg = nim_pkg,
          lambda = param_samples[sample_ind, c('lambda[1]', 'lambda[2]')], 
          betas_tx = betas_tx, 
          template_bins = template_bins,
          depth_threshold = deep_dive_depth,
          deeper = TRUE
        )
        
        # sampled time to deep depth  
        length(fwd_sim$depths) - 1
        
      })
      
      list(list(
        covariate_combination = covariate_combinations,
        samples = post_pred_samples
      ))
      
    }, 
    pattern = map(covariate_combinations), 
    deployment = 'worker',
    memory = 'transient'
  )
)
