validation_nstep_targets = list(
  
  # define "n" for the n-step ahead prediction distribution
  tar_target(validation_nstep, 2),
  
  # sample indices for validating the n-step ahead prediction distribution
  tar_target(
    name = validation_data_nstep,
    command = {
      
      nstep = validation_nstep
      
      # define proportion of data used for validation
      propval = .1
      
      # load source data
      nim_pkg = data_pkg
      
      # transform covariates for each segment
      nim_pkg$data$covariates = do.call(
        cbind, lapply(1:nrow(data_pkg$consts$segments), function(seg_ind) {
          seg = data_pkg$consts$segments[seg_ind,]
          covariate_tx(
            covariates = data_pkg$data$covariates[
              , seg['start_ind']:seg['end_ind']
            ],
            control = covariate_tx_control
          )
      }))
      
      # trim segment definitions to exclude timepoints with missing covariates
      # (i.e., when there is not enough history to compute lagged cov's.),
      # and also timepoints without enough forward observations for validation
      nim_pkg$consts$segments = t(
        apply(nim_pkg$consts$segments, 1, function(r) {
          # original segment indices
          seg_inds = r['start_ind']:r['end_ind']
          # determine which indices have fully-defined covariates
          defined_inds = complete.cases(t(nim_pkg$data$covariates[,seg_inds]))
          # subset segment indices to ensure covariates are always defined
          seg_inds = seg_inds[defined_inds]
          # subset segment indices to ensure n-step validation can be done
          seg_inds = seg_inds[1:(length(seg_inds)-nstep)]
          # update segment definition
          r['start_ind'] = min(seg_inds)
          r['end_ind'] = max(seg_inds)
          r['length'] = length(seg_inds)
          r
        })
      )
      
      # fully enumerate the points that can be used for validation
      feasible_val_points = do.call(c, 
        apply(nim_pkg$consts$segments, 1, function(r) {
          r['start_ind']:r['end_ind']
      }))
      
      # sample points to use for validation
      val_inds = sample(
        x = feasible_val_points, 
        size = length(feasible_val_points) * propval, 
        replace = FALSE
      )
      
      val_inds
    }
  ),
  
  # define validation batch job targets, to be run externally
  tar_target(
    name = validation_config_nstep, 
    command = {
      
      # number of sub-tasks in array job
      ntasks = 100
      
      #
      # load model outputs
      #
      
      # # load posterior samples and model package
      # tar_load(fit_marginalized_model)
      
      fit_marginalized_model = list(
        samples = "output/mcmc/fit_marginalized_model/samples.rds",
        package = "output/mcmc/fit_marginalized_model/nim_pkg.rds"
      )
      
      samples = readRDS(fit_marginalized_model$samples)

      # set burn-in
      burn = 1:(nrow(samples)*.1)
      
      # identify posterior samples that will be used in the posterior analysis
      posterior_sample_inds = (1:nrow(samples))[-burn]
      
      # enumerate all simulations to run
      validation_step_tgts = expand.grid(
        posterior_sample_ind = posterior_sample_inds,
        validation_nstep_ind = validation_data_nstep
      )
      
      # divide simulation requirements across number of array jobs
      task_groups_raw = parallel::splitIndices(
        nx = nrow(validation_step_tgts), ncl = ntasks
      )
      
      # steps that each task should process
      task_groups = lapply(task_groups_raw, function(inds) {
        validation_step_tgts[inds,]
      })
      
      # save outputs
      f = file.path('Rscripts', 'validation_nstep')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '_task_groups.rds', sep = ''))
      
      saveRDS(task_groups, f)
      
      f
    }
  )
   
)