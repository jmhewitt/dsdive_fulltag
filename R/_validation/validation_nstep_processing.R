validation_nstep_processing_script = tar_target(
  name = validation_nstep_processing, 
  command = {

    tar_load(fit_marginalized_model)
    
    # fit_marginalized_model = list(
    #   samples = "output/mcmc/fit_marginalized_model/samples.rds",
    #   package = "output/mcmc/fit_marginalized_model/nim_pkg.rds"
    # )
    
    # load model configuration
    nim_pkg = readRDS(fit_marginalized_model$package)
    
    # scrape output files
    f = dir(path = file.path('Rscripts', 'validation_nstep', 'samples'), 
            pattern = 'rds', full.names = TRUE)
    
    # verify we have the correct amount of data
    if(length(f) != n_validation_tasks) {
      stop('Number of validation files does not match validation task count.')
    }
    
    # load and combine all posterior predictive samples for validation
    samples = do.call(rbind, lapply(f, readRDS))
    
    # compute predictive distributions, validation info, and scores
    validation_ind = unique(samples$validation_nstep_ind)
    validation_distns = lapply(validation_ind, function(ind) {
      # only work with posterior samples for this validation point
      df = samples %>% filter(validation_nstep_ind == ind)
      # summarize predictive distribution and data
      res = list(
        validation_nstep_ind = ind,
        latent_state_training_inds = df$latent_state_training_inds[1],
        from_depth_bin = nim_pkg$data$depth_bins[ind],
        to_depth_bin_1 = nim_pkg$data$depth_bins[ind + 1],
        to_depth_bin_2 = nim_pkg$data$depth_bins[ind + 2],
        pred_depth_bin_1 = tabulate(
          bin = df$pred_depth_bin_1, 
          nbins = nim_pkg$consts$n_bins
        ) / nrow(df),
        pred_depth_bin_2 = tabulate(
          bin = df$pred_depth_bin_2, 
          nbins = nim_pkg$consts$n_bins
        ) / nrow(df)
      )
      # compute validation scores
      res$scores = data.frame(
        validation_nstep_ind = res$validation_nstep_ind,
        latent_state_training_inds = res$latent_state_training_inds,
        crps_bin_1 = crps_sample(
          y = res$to_depth_bin_1, dat = df$pred_depth_bin_1
        ),
        crps_bin_2 = crps_sample(
          y = res$to_depth_bin_2, dat = df$pred_depth_bin_2
        ),
        mae_bin_1 = abs(res$to_depth_bin_1 - which.max(res$pred_depth_bin_1)),
        mae_bin_2 = abs(res$to_depth_bin_2 - which.max(res$pred_depth_bin_2))
      )
      # return results
      res
    })
    
    validation_distns
  }
)