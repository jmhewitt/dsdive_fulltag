fit_marginalized_model_script = tar_target(
  name = fit_marginalized_model,
  command = {
    
    # we must mainly copy data_pkg contents into a local variable so that the 
    # targets package will recognize data_pkg as an external dependency, but 
    # it is also reasonably good practice to avoid modifying the data_pkg input
    nim_pkg = data_pkg
    
    #
    # transform and compress covariates
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

    # get unique combinations of covariates that drive movement
    covs = t(unique(t(do.call(cbind, apply(
      nim_pkg$consts$segments, 1, function(r) {
        nim_pkg$data$covariates[,r['start_ind']:r['end_ind']]
    })))))
    
    # compress covariates
    covariate_map = ht()
    for(i in 1:ncol(covs)) {
      tmp = covs[, i, drop = FALSE]
      rownames(tmp) = NULL
      covariate_map[tmp] = i
    }
    
    # add compressed id's to data
    nim_pkg$data$covariateId = sapply(
      1:ncol(nim_pkg$data$covariates), function(i) {
        tmp = nim_pkg$data$covariates[, i, drop = FALSE]
        rownames(tmp) = NULL
        id = covariate_map[tmp]
        ifelse(is.null(id), NA, id)
    })
    
    # add unique covariates and their count to data package
    nim_pkg$consts$covariates_unique = covs
    nim_pkg$consts$n_covariate_combinations = ncol(covs)
    nim_pkg$consts$n_covariates = nrow(covs)

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
    
    # treat infinitesimal depth bin transition probabilities as fixed;
    # values taken from Hewitt et al. (2021)
    nim_pkg$consts$pi = c(.97, .5, .05)
    
    # initial speed values and lambda dimension
    nim_pkg$inits$lambda = c(.3, 1.5)
    nim_pkg$consts$n_lambda_class = length(nim_pkg$inits$lambda)
    
    # initialize depth bin transition matrices
    nim_pkg$inits$depth_tx_mat = array(
      dim = c(nim_pkg$consts$n_stages, 
              nim_pkg$consts$n_bins, 
              nim_pkg$consts$n_bins)
    )
    for(i in 1:nim_pkg$consts$n_stages) {
      nim_pkg$inits$depth_tx_mat[i, , ] <- expm::expm(
        nim_pkg$consts$tstep * buildInfinitesimalGenerator(
          pi = nim_pkg$consts$pi[nim_pkg$consts$stage_defs[i, 1]],
          lambda = nim_pkg$inits$lambda[nim_pkg$consts$stage_defs[i, 2]],
          M = nim_pkg$consts$n_bins,
          stage = 6,
          widths = nim_pkg$consts$widths
        )
      )
    }
    
    # initialize stage transition matrices
    nim_pkg$inits$stage_tx_mat = array(
      data = 1/nim_pkg$consts$n_stages,
      dim = c(nim_pkg$consts$n_subjects,
              nim_pkg$consts$n_stages,
              nim_pkg$consts$n_stages,
              nim_pkg$consts$n_covariate_combinations)
    )
    
    # initial stage transition covariate effects per individual (random effects)
    nim_pkg$inits$beta_tx = array(
      data = 0,
      dim = c(nim_pkg$consts$n_subjects,
              nim_pkg$consts$n_covariates,
              nim_pkg$consts$n_stages,
              nim_pkg$consts$n_stages - 1)
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
    
    # hyperprior for population-level effect centers
    nim_pkg$consts$beta_tx_mu_prior_mean = numeric(
      length = nim_pkg$consts$n_stages - 1
    )
    nim_pkg$consts$beta_tx_mu_prior_cov = diag(
      x = 1e2, 
      nrow = nim_pkg$consts$n_stages - 1, 
      ncol = nim_pkg$consts$n_stages - 1
    )
    
    # hyperprior for population-level effect precisions
    nim_pkg$consts$beta_tx_prec_prior = prec_mat
    nim_pkg$consts$beta_tx_prec_prior_k = nim_pkg$consts$n_stages
    
    #
    # nimble model fitting
    #
    
    # TRUE to estimate separate covariate effects for each individual
    nim_pkg$consts$random_effects = TRUE
    
    # remove objects not used in model
    # nim_pkg$data$covariates = NULL
    nim_pkg$data$times = NULL
    
    mod = nimbleModel(
      code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data,
      inits = nim_pkg$inits
    )
    
    cmod = compileNimble(mod)
    
    # verify model has a finite likelihood
    if(!is.finite(cmod$calculate())) {
      stop('Initial likelihood is not finite')
    }
    
    conf = configureMCMC(mod)
    
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
    
    conf$addMonitors2(c('beta_tx', 'depth_tx_mat'))
    
    mcmc = buildMCMC(conf)
    
    cmcmc = compileNimble(mcmc)
    
    niter = 1e4
    
    fpath = file.path('output', 'mcmc', tar_name())
    dir.create(path = fpath, showWarnings = FALSE, recursive = TRUE)
    
    # save model package
    f_pkg = file.path(fpath, 'nim_pkg.rds')
    saveRDS(nim_pkg, file = f_pkg)
    
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
    
    
    
    # some exploration of model output
    # 
    # burn = 1:1e3
    # plot(samples[,'beta_tx[5, 5, 4]'])
    # plot(mcmc(samples[-burn,'lambda[1]']))
    # plot(mcmc(samples[-burn,'lambda[2]']))
    # sort(effectiveSize(mcmc(samples[-burn,])))
    # 
    # 1-rejectionRate(mcmc(samples[-burn,'lambda[1]']))
    
    list(
      samples = fpath,
      package = f_pkg
    )
  }
)