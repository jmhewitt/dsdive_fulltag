fit_model = function(data_pkg, covariate_tx_control, movement_classes,
                     random_inits, out_dir, random_betas) {
  # Parameters:
  #   random_inits - TRUE to sample initial model parameters from prior
  #   random_betas - FALSE to start all betas from 0 and cov mat's from I
  #   out_dir - directory in which to save posterior samples and configuration
  
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
  # nimble model fitting
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
  
  # sample initial model parameter values from the prior
  if(random_inits) {
    initializing_params = TRUE
    while(initializing_params) { 
      # resampling overcomes the possibility of choosing bad parameter values,
      # in particular for lambda that are both approximately 0
      message('(Re-)Sampling initial parameters')
      # draw top-level parameters
      cmod$simulate(nodes = c('pi','lambda'))
      if(random_betas) {
        cmod$simulate(nodes = c('beta_tx_mu', 'beta_tx_prec'))
      }
      # recall pi[2] is fixed
      cmod$pi[2] = .5
      # ensure lambda constraints are met
      cmod$lambda = sort(cmod$lambda)
      # update dependent model parameters and components
      cmod$calculate()
      # draw random effects
      if(random_betas) {
        cmod$simulate(nodes = 'beta_tx')
      }
      # update dependent model parameters and components
      initializing_params = !is.finite(cmod$calculate())
    }
  }
  
  # verify model has a finite likelihood
  if(!is.finite(cmod$calculate())) {
    stop('Initial likelihood is not finite')
    # debugging support to find the nodes with poor definition
    ll = sapply(cmod$getNodeNames(), function(nom) {
      cmod$calculate(nom)
    })
    names(ll[which(!is.finite(ll))])
  }
  
  message('Building MCMC sampler')
  
  conf = configureMCMC(mod)
  
  # treat directional preference for free movement as known
  conf$removeSampler('pi[2]')
  
  # NOTE: The alternate samplers seem to slow down sampling rates by 10x while
  #  at most doubling the average effective sample size, and not changing the 
  #  worst effective sampling size substantially.
  # 
  # # get top-level names and groupings of variables being sampled
  # sampling_targets = sapply(conf$getSamplers(), function(x) x$target)
  # sampling_target_groups = gsub(
  #   pattern = '\\[.*\\]', 
  #   replacement = '', 
  #   x = sampling_targets
  # )
  # sampling_groups = unique(sampling_target_groups)
  # 
  # # replace block samplers over gaussian nodes with elliptical slice samplers
  # for(tgt_group in c('beta_tx_mu', 'beta_tx')) {
  #   if(tgt_group %in% sampling_groups) {
  #     conf$removeSamplers(tgt_group)
  #     for(tgt in sampling_targets[sampling_target_groups == tgt_group]) {
  #       conf$addSampler(target = tgt, type = 'ess')
  #     }
  #   }
  # }
  # 
  # # replace univariate RW samplers with slice samplers
  # for(tgt_group in c('lambda')) {
  #   if(tgt_group %in% sampling_groups) {
  #     conf$removeSamplers(tgt_group)
  #     for(tgt in sampling_targets[sampling_target_groups == tgt_group]) {
  #       conf$addSampler(target = tgt, type = 'slice')
  #     }
  #   }
  # }
  
  # save derived quantities that are required for posterior sampling
  conf$addMonitors2(c('beta_tx', 'depth_tx_mat'))
  
  mcmc = buildMCMC(conf)
  
  cmcmc = compileNimble(mcmc)
  
  niter = 1e4
  
  fpath = out_dir
  
  # make sure output directory is empty (i.e., clear previous model output)
  unlink(x = fpath, recursive = TRUE)
  
  # (re-)create output directory
  dir.create(path = fpath, showWarnings = FALSE, recursive = TRUE)
  
  # save model package
  f_pkg = file.path(fpath, 'nim_pkg.rds')
  saveRDS(nim_pkg, file = f_pkg)
  
  message('Running MCMC sampler')
  
  # run model, saving output along the way
  runCheckpointMCMC(
    mcmc = cmcmc, 
    nsamples = niter, 
    out.path = fpath,
    out.name = 'samples', 
    out.size = 1024, 
    out.time = 3600/2,
    verbose = TRUE
  )
  
  list(
    samples = fpath,
    package = f_pkg
  )
}
