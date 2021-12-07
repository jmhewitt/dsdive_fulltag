library(targets)
library(nimble)
library(coda)

# load custom likelihood
source(
  dir(path = file.path('R', 'util', 'statespace_tools'), 
      pattern = '.R', full.names = TRUE)
)

# load matrix exponentials
source(
  dir(path = file.path('R', 'util', 'expokit'), 
      pattern = '.R', full.names = TRUE)
)

# load nimble components
lapply(
  dir(path = file.path('R', 'nimble'), pattern = '.R', full.names = TRUE),
  source
)

# load hashmap function
source(file.path('R', 'util', 'hashmap.R'))

tar_load(nim_pkg)

# get unique combinations of covariates that drive movement
covs = t(unique(t(nim_pkg$data$covariates)))

# compress covariates
covariate_map = ht()
for(i in 1:ncol(covs)) {
  tmp = covs[, i, drop = FALSE]
  rownames(tmp) = NULL
  covariate_map[tmp] = i
}
nim_pkg$data$covariateId = sapply(1:ncol(nim_pkg$data$covariates), function(i) {
  tmp = nim_pkg$data$covariates[, i, drop = FALSE]
  rownames(tmp) = NULL
  covariate_map[tmp]
})

#
# add additional package elements
#


nim_pkg$consts$covariates_unique = covs

nim_pkg$consts$init_stages = rep(
  1/nim_pkg$consts$n_stages, nim_pkg$consts$n_stages
)

nim_pkg$consts$n_covariate_combinations = ncol(covs)

nim_pkg$consts$pi = c(.01, .5, .99)
nim_pkg$inits$lambda = rep(.5, nim_pkg$consts$n_lambda_class)
nim_pkg$consts$tstep = 300
nim_pkg$consts$widths = tar_read(template_bins)$halfwidth * 2

# initialize depth bin transition matrices
nim_pkg$inits$depth_tx_mat = array(
  dim = c(nim_pkg$consts$n_stages, nim_pkg$consts$n_bins, nim_pkg$consts$n_bins)
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
  dim = c(nim_pkg$consts$n_stages, nim_pkg$consts$n_stages, 
          nim_pkg$consts$n_covariate_combinations)
)
for(i in 1:nim_pkg$consts$n_covariate_combinations) {
  nim_pkg$inits$stage_tx_mat[, , i] <- stageTxMats(
    betas = nim_pkg$inits$beta_tx,
    covariates = nim_pkg$consts$covariates_unique[, i, drop = FALSE],
    n_timepoints = 1
  )[, , 1]
}


#
# nimble modeling
#


buildBinTxMats = nimbleFunction(
  run = function(n_stages = integer(0),
                 pi = double(1),
                 lambda = double(1),
                 n_bins = integer(0),
                 widths = double(1),
                 tstep = double(0),
                 stage_defs = double(2)) {
    
    returnType(double(3))
    
    depth_tx_mat <- array(dim = c(n_stages, n_bins, n_bins), init = FALSE)
    H <- matrix(nrow = n_bins, ncol = n_bins, init = FALSE)
    
    for(i in 1:n_stages) {
      depth_tx_mat[i, 1:n_bins, 1:n_bins] <- expocall_gpadm(
        H = buildInfinitesimalGenerator(
          pi = pi[stage_defs[i, 1]],
          lambda = lambda[stage_defs[i, 2]],
          M = n_bins,
          stage = 6,
          widths = widths[1:n_bins]
        ),
        t = tstep,
        nrows = n_bins,
        ncols = n_bins
      )[1:n_bins, 1:n_bins]
    }
    
    return(depth_tx_mat)
  }
)


# model!
modelCode = nimble::nimbleCode({
  
  # speed classes
  for(j in 1:n_lambda_class) {
    lambda[j] ~ dgamma(shape = .01, rate = .01)
  }
  
  # constraint to assist in identifying speed classes
  constraint_data ~ dconstraint(lambda[1] <= lambda[2])
  
  # stage transition matrix effects
  for(i in 1:n_covariates) {
    for(j in 1:n_stages) { # stage being transitioned from
      beta_tx[i,j,1:(n_stages-1)] ~ dmnorm(mean = beta_tx_mu[1:(n_stages-1)],
                                           cov = beta_tx_cov[1:(n_stages-1),
                                                             1:(n_stages-1)])
    }
  }
  
  # stage transition matrices
  for(i in 1:n_covariate_combinations) {
    stage_tx_mat[1:n_stages, 1:n_stages, i] <- stageTxMats(
      betas = beta_tx[1:n_covariates, 1:n_stages, 1:(n_stages-1)],
      covariates = covariates_unique[1:n_covariates, i:i],
      n_timepoints = 1
    )[1:n_stages, 1:n_stages, 1]
  }
  
  # depth bin transition matrices
  for(i in 1:n_stages) {
    depth_tx_mat[i, 1:n_bins, 1:n_bins] <- expocall_gpadm(
      H = buildInfinitesimalGenerator(
        pi = pi[stage_defs[i, 1]],
        lambda = lambda[stage_defs[i, 2]],
        M = n_bins,
        stage = 6,
        widths = widths[1:n_bins]
      ),
      t = tstep,
      nrows = n_bins,
      ncols = n_bins
    )[1:n_bins, 1:n_bins]
  }
  
  # largest sampling unit is a sequence of depth bins
  for(seg_num in 1:n_segments) {
    # likelihood for depth bin observations
    depths[segments[seg_num,1]:segments[seg_num,4]] ~ dstatespace(
      obs_lik_dict = depth_tx_mat[1:n_stages, 1:n_bins, 1:n_bins],
      txmat_dict = stage_tx_mat[
        1:n_stages, 1:n_stages, 1:n_covariate_combinations
      ],
      txmat_seq = covariateId[segments[seg_num,1]:segments[seg_num,4]],
      x0 = init_stages[1:n_stages],
      num_obs_states = n_bins,
      num_latent_states = n_stages,
      nt = segments[seg_num,2]
    )
  }
  
})

mod = nimbleModel(
  code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data, 
  inits = nim_pkg$inits
)

cmod = compileNimble(mod)

cmod$calculate()

conf = configureMCMC(mod)

mcmc = buildMCMC(conf)

cmcmc = compileNimble(mcmc)

niter = 1e3
tick = proc.time()[3]
cmcmc$run(niter = niter)
tock = proc.time()[3]
niter / (tock - tick)

burn = 1:700
samples = as.matrix(cmcmc$mvSamples)
plot(samples[,'beta_tx[5, 5, 4]'])
plot(samples[,'lambda[1]'])
plot(mcmc(samples[-burn,'lambda[2]']))
sort(effectiveSize(mcmc(samples[-burn,])))

1-rejectionRate(mcmc(samples[-burn,'beta_tx[5, 5, 4]']))
