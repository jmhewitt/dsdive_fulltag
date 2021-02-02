fit = function(nim_pkg) {
  browser()
  
  library(targets)
  
  source('_targets.R')
  
  
  tar_load(nim_pkg)
  tar_load(template_bins)
  
  alpha = rbind(
    qlogis(c(
      .99, # deep_descent
      .5,  # deep_forage
      .01, # deep_ascent
      .97, # shallow_descent
      .03, # shallow_ascent
      .5   # free_surface
    )),
    rep(0, 6), # slopes to ignore depth in common stage/bin tx covariate matrix
    rep(0, 6), # slopes to ignore deep depth indicator
    rep(0, 6) # slopes to ignore shallow depth indicator
  )
  
  # intercepts for bin transition vertical speeds
  beta = rbind(
    log(c(
      1.6,  # deep_descent
      .03,  # deep_forage
      1.4,  # deep_ascent
      .6,   # shallow_descent
      .6,   # shallow_ascent
      .2    # free_surface
    )),
    rep(0, 6), # slopes to ignore depth in common stage/bin tx covariate matrix
    rep(0, 6), # slopes to ignore deep depth indicator
    rep(0, 6) # slopes to ignore shallow depth indicator
  )
  
  # covariates
  covariates = rbind(
    intercept = rep(1, length(nim_pkg$data$depths)),
    depth = template_bins$center[nim_pkg$data$depths],
    deep_depth = template_bins$center[nim_pkg$data$depths] >= 800,
    shallow_depth = template_bins$center[nim_pkg$data$depths] < 800 
  )
  
  # stage transition matrix, from which we will derive multinomial logit coefs.
  o = c('deep_descent', 'deep_forage', 'deep_ascent', 'shallow_descent', 
        'shallow_ascent', 'free_surface')
  ideal_stage_txmat = matrix(0, nrow = 6, ncol = 6, 
                             dimnames = list(rownames = o, colnames = o))
  ideal_stage_txmat['deep_descent', 'deep_forage'] = 1/15
  ideal_stage_txmat['deep_forage', 'deep_ascent'] = 1/20
  ideal_stage_txmat['deep_ascent', 'deep_descent'] = 1/15 * 1/5 * 95/100
  ideal_stage_txmat['deep_ascent', 'shallow_descent'] = 1/15 * 4/5
  ideal_stage_txmat['deep_ascent', 'free_surface'] = 1/15 * 1/5 * 5/100
  ideal_stage_txmat['shallow_descent', 'shallow_ascent'] = 1/7
  ideal_stage_txmat['shallow_ascent', 'deep_descent'] = 1/7 * 1/5 * 95/100
  ideal_stage_txmat['shallow_ascent', 'shallow_descent'] = 1/7 * 4/5
  ideal_stage_txmat['shallow_ascent', 'free_surface'] = 1/7 * 1/5 * 5/100
  ideal_stage_txmat['free_surface', 'deep_descent'] = 4/5 * 9/10 # 1/15 * 4/5
  ideal_stage_txmat['free_surface', 'shallow_descent'] = 1/5 * 9/10 # 1/15 * 1/5
  # ideal_stage_txmat = ideal_stage_txmat / 100
  diag(ideal_stage_txmat) = 1 - rowSums(ideal_stage_txmat)
  
  # multinomial logit coefs associated with ideal stage transition matrix
  betas_tx = rbind(
    intercept = c(
      log(ideal_stage_txmat['deep_descent', 'deep_forage']) - 
        log(ideal_stage_txmat['deep_descent', 'deep_descent']),
      log(ideal_stage_txmat['deep_forage', 'deep_ascent']) - 
        log(ideal_stage_txmat['deep_forage', 'deep_forage']),
      log(ideal_stage_txmat['deep_ascent', 'deep_descent']) - 
        log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
      log(ideal_stage_txmat['deep_ascent', 'shallow_descent']) - 
        log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
      log(ideal_stage_txmat['deep_ascent', 'free_surface']) - 
        log(ideal_stage_txmat['deep_ascent', 'deep_ascent']),
      log(ideal_stage_txmat['shallow_descent', 'shallow_ascent']) - 
        log(ideal_stage_txmat['shallow_descent', 'shallow_descent']),
      log(ideal_stage_txmat['shallow_ascent', 'deep_descent']) - 
        log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
      log(ideal_stage_txmat['shallow_ascent', 'shallow_descent']) - 
        log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
      log(ideal_stage_txmat['shallow_ascent', 'free_surface']) - 
        log(ideal_stage_txmat['shallow_ascent', 'shallow_ascent']),
      log(ideal_stage_txmat['free_surface', 'deep_descent']) - 
        log(ideal_stage_txmat['free_surface', 'free_surface']),
      log(ideal_stage_txmat['free_surface', 'shallow_descent']) - 
        log(ideal_stage_txmat['free_surface', 'free_surface'])
    ),
    depth = rep(0, 11),
    deep_depth =  0 * c( # Use -Inf to disallow certain stage tx's. at deep depths
      0,    # deep_descent -> deep_forage
      0,    # deep_forage -> deep_ascent 
      0, # deep_ascent -> deep_descent
      -1e2, # deep_ascent -> shallow_descent
      -1e2,    # deep_ascent -> free_surface
      1e2, # shallow_descent -> shallow_ascent
      1e2, # shallow_ascent -> deep_descent
      -1e2, # shallow_ascent -> shallow_descent
      -1e2, # shallow_ascent -> free_surface
      1e2,    # free_surface -> deep_descent
      1e2  # free_surface -> shallow_descent; at-depth descent is always deep_*
    ),
    shallow_depth = 0 * c( # Use -Inf to disallow certain stage tx's. at deep depths
      -1e2,    # deep_descent -> deep_forage
      1e2,    # deep_forage -> deep_ascent 
      0, # deep_ascent -> deep_descent
      0,    # deep_ascent -> shallow_descent; possible if obs. is at end of dive
      0,# deep_ascent -> free_surface
      0,  # shallow_descent -> shallow_ascent
      0,    # shallow_ascent -> deep_descent
      0, # shallow_ascent -> shallow_descent
      0, # shallow_ascent -> free_surface
      0,    # free_surface -> deep_descent
      0  # free_surface -> shallow_descent; at-depth descent is always deep_*
    )
  )
  
  nim_pkg$inits = list(
    # population-level effects
    alpha_mu = alpha,
    beta_mu = beta,
    alpha_var = matrix(1, nrow = nrow(alpha), ncol = ncol(alpha)),
    beta_var = matrix(1, nrow = nrow(beta), ncol = ncol(beta)),
    # individual-level random effects
    alpha = array(
      data = alpha, 
      dim = c(nrow(alpha), ncol(alpha), nim_pkg$consts$n_subjects)
    ),
    beta = array(
      data = beta, 
      dim = c(nrow(beta), ncol(beta), nim_pkg$consts$n_subjects)
    ),
    betas_tx = array(
      data = betas_tx, 
      dim = c(nrow(betas_tx), ncol(betas_tx), nim_pkg$consts$n_subjects)
    )
  )
  
  nim_pkg$data$covariates = covariates
  
  nim_pkg$consts$n_covariates = nrow(nim_pkg$data$covariates)
  
  nim_pkg$consts$n_txmat_entries = length(nim_pkg$consts$transition_matrices)
  nim_pkg$consts$n_txmat_types = nrow(nim_pkg$consts$pi_discretization)
  
  nim_pkg$consts$n_stages = length(nim_pkg$consts$movement_types)
  
  nim_pkg$consts$n_stage_txs = ncol(betas_tx)
  
  nim_pkg$consts$n_timepoints = length(nim_pkg$data$depths)
  
  nim_pkg$consts$segments = cbind(
    nim_pkg$consts$segments,
    end_ind = nim_pkg$consts$segments[, 'start_ind'] + 
      nim_pkg$consts$segments[, 'length'] - 
      1
  ) 
  
  nim_pkg$data$transition_matrices = nim_pkg$consts$transition_matrices
  nim_pkg$consts$transition_matrices = NULL
  
  source('R/nimble/model.R')
  
  model = nimbleModel(code = modelCode, constants = nim_pkg$consts, 
                      data = nim_pkg$data, inits = nim_pkg$inits, 
                      name = 'fulltag')
  
  model_c = compileNimble(model, projectName = 'fulltag')
  
  conf = configureMCMC(model)
  
  # remove samplers for movement parameter effects not being used
  for(i in 2:nim_pkg$consts$n_covariates) {
    for(j in 1:nim_pkg$consts$n_stages) {
      conf$removeSampler(paste('alpha_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('beta_var[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('beta_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('alpha_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        conf$removeSampler(paste('alpha[', i, ', ', j, ', ', k, ']', sep = ''))
        conf$removeSampler(paste('beta[', i, ', ', j, ', ', k, ']', sep = ''))
      }
    }
  }
  
  # only the first row
  conf$removeSampler('betas_tx')
  

  # latent stage vector samplers
  conf$removeSampler('stages')
  for(seg_ind in 1:nim_pkg$consts$n_segments) {
    start_ind = nim_pkg$consts$segments[seg_ind, 'start_ind']
    end_ind = start_ind + nim_pkg$consts$segments[seg_ind, 'length'] - 1
    conf$addSampler(
      target = paste('stages[', start_ind, ':', end_ind, ']', sep = ''),
      type = 'Stage',
      control = list(
        betas_tx_node = paste('betas_tx'),
        stage_supports = nim_pkg$data$stage_supports[, start_ind:end_ind]
      )
    )
  }
  
  conf$addMonitors('stages')
  
  mcmc = buildMCMC(conf)
  
  mcmc_c = compileNimble(mcmc, resetFunctions = TRUE, showCompilerOutput = TRUE,
                         projectName = 'fulltag')
  
  mcmc_c$run(niter = 10)
  samples = as.matrix(mcmc_c$mvSamples)
  pryr::object_size(samples)
  model$calculate()
  
  str(as.matrix(mcmc$mvSamples))
  
}