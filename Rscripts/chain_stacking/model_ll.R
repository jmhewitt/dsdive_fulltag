model_ll = function(data_pkg, covariate_tx_control, movement_classes,
                    samples_dir, out_dir) {
  # Parameters:
  #   samples_dir - directory from which to load  posterior samples and config.
  #   out_dir - directory in which to save likelihood computations
  
  # we build the model config/data from the input
  nim_pkg = data_pkg
  
  #
  # transform covariates
  #
  
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
  # (i.e., as occur when there is not enough history to compute lagged cov's.)
  nim_pkg$consts$segments = t(
    apply(nim_pkg$consts$segments, 1, function(r) {
      # original segment indices
      seg_inds = r['start_ind']:r['end_ind']
      # determine which indices have fully-defined covariates
      defined_inds = complete.cases(t(nim_pkg$data$covariates[,seg_inds]))
      # subset segment indices
      seg_inds = seg_inds[defined_inds]
      # update segment definition
      r['start_ind'] = min(seg_inds)
      r['end_ind'] = max(seg_inds)
      r['length'] = length(seg_inds)
      r
    })
  )
  
  # add covariate count to package
  nim_pkg$consts$n_covariates = nrow(nim_pkg$data$covariates)
  
  #
  # add additional nimble package elements
  #
  
  # storage for initial values of nuisance parameters, etc.
  nim_pkg$inits = list()
  
  # support ordering estimated speed parameters
  nim_pkg$data$constraint_data = 1
  
  # movement stages to model
  nim_pkg$consts$stage_defs = movement_classes$stage_defs
  nim_pkg$consts$n_stages = nrow(nim_pkg$consts$stage_defs)
  
  # marginal likelihood for each segment uses this prior as the initial 
  # distribution of latent movement classes
  nim_pkg$consts$init_stages = rep(
    1/nim_pkg$consts$n_stages, nim_pkg$consts$n_stages
  )
  
  # constrain support for pi's to help identify the stages as ascent/descent
  nim_pkg$consts$pi_prior_min = c(.5, 0, 0)
  nim_pkg$consts$pi_prior_max = c(1, 1, .5)
  
  # initial values for infinitesimal depth bin transition probabilities taken
  # from Hewitt et al. (2021); nim_pkg$inits$pi[2] will be treated as fixed
  nim_pkg$inits$pi = c(.97, .5, .05)
  nim_pkg$consts$n_pi_class = length(nim_pkg$inits$pi)
  
  # initial speed values and lambda dimension
  nim_pkg$consts$n_lambda_class = length(unique(
    nim_pkg$consts$stage_defs[,'lambda_ind']
  ))
  nim_pkg$inits$lambda = c(.3, 1.5)
  
  # hyperprior for population-level effect centers
  nim_pkg$consts$beta_tx_mu_prior_mean = numeric(
    length = nim_pkg$consts$n_stages - 1
  )
  nim_pkg$consts$beta_tx_mu_prior_cov = diag(
    x = 1e2, 
    nrow = nim_pkg$consts$n_stages - 1, 
    ncol = nim_pkg$consts$n_stages - 1
  )
  
  # group level prior means or hierarchical centering for covariate effects
  nim_pkg$inits$beta_tx_mu = array(
    data = 0,
    dim = c(nim_pkg$consts$n_covariates,
            nim_pkg$consts$n_stages,
            nim_pkg$consts$n_stages - 1)
  )
  
  # group level prior precision matrices for covariate effects
  prec_mat = diag(
    x = 1e-2, 
    nrow = nim_pkg$consts$n_stages - 1,
    ncol = nim_pkg$consts$n_stages - 1
  )
  nim_pkg$inits$beta_tx_prec = array(
    dim = c(nim_pkg$consts$n_covariates,
            nim_pkg$consts$n_stages,
            nim_pkg$consts$n_stages - 1,
            nim_pkg$consts$n_stages - 1)
  )
  for(i in 1:nim_pkg$consts$n_covariates) {
    for(j in 1:nim_pkg$consts$n_stages) {
      nim_pkg$inits$beta_tx_prec[i,j,,] = prec_mat
    }
  }
  
  # hyperprior for population-level effect precisions
  nim_pkg$consts$beta_tx_prec_prior = prec_mat
  nim_pkg$consts$beta_tx_prec_prior_k = nim_pkg$consts$n_stages
  
  # initial stage transition covariate effects per individual (random effects)
  nim_pkg$inits$beta_tx = array(
    data = 0,
    dim = c(nim_pkg$consts$n_subjects,
            nim_pkg$consts$n_covariates,
            nim_pkg$consts$n_stages,
            nim_pkg$consts$n_stages - 1)
  )
  
  # initialize depth bin transition matrices
  nim_pkg$inits$depth_tx_mat = array(
    dim = c(nim_pkg$consts$n_stages, 
            nim_pkg$consts$n_bins, 
            nim_pkg$consts$n_bins)
  )
  for(i in 1:nim_pkg$consts$n_stages) {
    nim_pkg$inits$depth_tx_mat[i, , ] <- expm::expm(
      nim_pkg$consts$tstep * buildInfinitesimalGenerator(
        pi = nim_pkg$inits$pi[nim_pkg$consts$stage_defs[i, 1]],
        lambda = nim_pkg$inits$lambda[nim_pkg$consts$stage_defs[i, 2]],
        M = nim_pkg$consts$n_bins,
        stage = 6,
        widths = nim_pkg$consts$widths
      )
    )
  }
  
  #
  # nimble model construction
  #
  
  # specify which covariates are modeled as random vs. fixed effects:
  #   - only model intercept terms as random effects
  nim_pkg$consts$ranef_inds = which(
    rownames(nim_pkg$data$covariates) %in% c('intercept')
  )
  nim_pkg$consts$fixef_inds = setdiff(
    1:nim_pkg$consts$n_covariates, nim_pkg$consts$ranef_inds
  )
  nim_pkg$consts$n_fixefs = length(nim_pkg$consts$fixef_inds)
  nim_pkg$consts$n_ranefs = length(nim_pkg$consts$ranef_inds)
  
  # ensure nimble interprets the covariate index vectors as vectors;
  # if n_fixefs or n_ranefs == 1, then nimble will try to interpret as scalar
  nim_pkg$consts$fixef_inds = c(nim_pkg$consts$fixef_inds, NA)
  nim_pkg$consts$ranef_inds = c(nim_pkg$consts$ranef_inds, NA)
  
  # TRUE to estimate separate/group covariate effects for each individual
  nim_pkg$consts$random_effects = nim_pkg$consts$n_ranefs > 0
  nim_pkg$consts$fixed_effects = nim_pkg$consts$n_fixefs > 0
  
  # remove objects not used in model
  # nim_pkg$data$covariates = NULL
  nim_pkg$data$times = NULL
  
  mod = nimbleModel(
    code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data, 
    inits = nim_pkg$inits, calculate = FALSE
  )
  
  cmod = compileNimble(mod)
  
  #
  # load, label, and merge posterior samples
  #

  mvSample_files = dir(
    path = samples_dir,
    pattern = 'mvSamples_[0-9]+',
    full.names = TRUE
  )

  mvSample2_files = dir(
    path = samples_dir,
    pattern = 'mvSamples2_[0-9]+',
    full.names = TRUE
  )

  samples = do.call(rbind, lapply(mvSample_files, readRDS))
  samples2 = do.call(rbind, lapply(mvSample2_files, readRDS))

  colnames(samples) = readRDS(dir(
    path = samples_dir,
    pattern = 'mvSamples_colnames',
    full.names = TRUE
  ))

  colnames(samples2) = readRDS(dir(
    path = samples_dir,
    pattern = 'mvSamples2_colnames',
    full.names = TRUE
  ))

  samples = cbind(samples, samples2)
  rm(samples2)
  
  # get top-level names and groupings of variables being sampled
  sampling_targets = colnames(samples)
  sampling_target_groups = gsub(
    pattern = '\\[.*\\]',
    replacement = '',
    x = sampling_targets
  )
  sampling_groups = unique(sampling_target_groups)
  
  #
  # evaluate likelihood
  #
  
  posterior_sample_inds = (1:nrow(samples))
  
  # identify model parameters (i.e., a "theta" vector)
  param_tgts = cmod$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  
  # identify model data (i.e., a "y_i" vector)
  data_tgts = data_tgts = grep(
    pattern = 'depth_bins', 
    x = cmod$getNodeNames(dataOnly = TRUE), 
    value = TRUE
  )
  
  ll = do.call(rbind, lapply(posterior_sample_inds, function(ind) {
    
    # transfer posterior sample of model parameters to model object
    for(tgt_group in sampling_groups) {
      cmod[[tgt_group]] = samples[
        ind, sampling_targets[sampling_target_groups == tgt_group]
      ]
    }
    
    # update model internals
    cmod$calculate()
    
    # extract likelihood components
    c(
      cmod$calculate(param_tgts),
      unname(sapply(data_tgts, function(tgt) cmod$calculate(tgt)))
    )
  }))
  
  colnames(ll) = c('theta', data_tgts)
  
  #
  # save outputs
  #
  
  fpath = out_dir
  
  # (re-)create output directory
  dir.create(path = fpath, showWarnings = FALSE, recursive = TRUE)
  
  fname = 'model_ll_samples.rds'
  
  saveRDS(ll, file = file.path(fpath, fname))
  
  0
}
