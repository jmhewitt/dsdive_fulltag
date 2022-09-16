#
# load workflow dependencies and packages
#

source('_targets.R')

sapply(tar_option_get('packages'), function(x) {
  library(x, character.only = TRUE)
})

tar_load(
  c('cape_hatteras_loc', 'covariate_tx_control', 'validation_nstep',
    'data_pkg', 'template_bins')
)


#
# load target number to process and associated data
#

taskGroup = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(is.na(taskGroup)) {
  taskGroup = 1
}

tasks = readRDS(
  file.path('Rscripts', 'validation_nstep', 
            'validation_config_nstep_task_groups.rds')
)[[taskGroup]]


#
# load model outputs
#

tar_load(fit_marginalized_model)

# fit_marginalized_model = list(
#   samples = "output/mcmc/fit_marginalized_model/samples.rds",
#   package = "output/mcmc/fit_marginalized_model/nim_pkg.rds"
# )

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
# re-build model so that we can use it to generate transition matrices
#

# uncompiled model
mod = nimbleModel(
  code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data,
  inits = nim_pkg$inits, calculate = FALSE
)

# compile model
cmod = compileNimble(mod)

# verify model has a finite likelihood
cmod$calculate()


#
# composition sample!
#

# identify the posterior samples to loop over
posterior_sample_inds = unique(tasks$posterior_sample_ind)

# draw samples for tasks by looping over the posterior samples to work with
validation_pred_samples = do.call(rbind, lapply(posterior_sample_inds, 
                                                function(post_ind) {
                                                  
  # transfer posterior sample of model parameters to model object
  for(tgt_group in sampling_groups) {
    cmod[[tgt_group]] = samples[
      post_ind, sampling_targets[sampling_target_groups == tgt_group]
    ]
  }
  
  # update model components
  cmod$calculate()
  
  # all of the validation timepoints that need to be processed for this sample
  validation_inds = tasks %>% 
    filter(posterior_sample_ind == post_ind) %>%
    select(validation_nstep_ind) %>% 
    unlist()
  
  # draw predictive samples
  pred_samples = do.call(rbind, lapply(validation_inds, function(val_ind) {

    # segment to which prediction time belongs
    seg_num = max(which(nim_pkg$consts$segments[,'start_ind'] <= val_ind))
    
    # indices used to initialize the latent state for predictions
    seg_inds = nim_pkg$consts$segments[seg_num,'start_ind']:val_ind
    
    # indices used to initialize the depths and covariates for predictions
    full_seg_inds = data_pkg$consts$segments[seg_num,'start_ind']:val_ind
    
    # extract subject id for segment
    subject_id = nim_pkg$consts$segments[seg_num, 'subject_id']
    
    # posterior predictive distribution for latent state before predictions
    xf = finalPred2LayerCompressed(
      obs_lik_dict = cmod$depth_tx_mat,
      obs = cmod$depth_bins[seg_inds] - 1,
      txmat_dict = cmod$stage_tx_mat[subject_id, , , ],
      txmat_seq = cmod$covariateId[seg_inds] - 1,
      x0 = cmod$init_stages,
      num_obs_states = nim_pkg$consts$n_bins
    ) 
    
    # simulate dive until a deep depth is reached
    sim = fwd_sim_dive(
      stage = sample(
        x = 1:nim_pkg$consts$n_stages, size = 1, prob = xf
      ),
      depth_bins = nim_pkg$data$depth_bins[full_seg_inds],
      covariates = data_pkg$data$covariates[, full_seg_inds],
      n_timepoints = validation_nstep,
      nim_pkg = nim_pkg,
      model = cmod,
      lon = cape_hatteras_loc['lon'],
      lat = cape_hatteras_loc['lat'],
      depth_threshold = Inf, # Don't stop sampling if depth thresh. exceeded
      template_bins = template_bins,
      subj = subject_id,
      covariate_tx_control = covariate_tx_control
    )
    
    # package results
    df = data.frame(
      validation_nstep_ind = val_ind,
      posterior_sample_ind = post_ind,
      latent_state_training_inds = length(seg_inds)
    )
    for(i in 1:validation_nstep) {
      df[[paste('pred_depth_bin_', i, sep ='')]] = sim$depth_bins[i]
    }
    
    df
  }))

}))


# clear extraneous rownames
rownames(validation_pred_samples) = NULL

# save output
f = file.path('Rscripts', 'validation_nstep', 'samples')
dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  validation_pred_samples, 
  file = file.path(f, paste('predictive_samples_', taskGroup, '.rds', sep = ''))
)
