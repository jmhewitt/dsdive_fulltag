sim_summary = function(sim_pkg, pi, sample_file, nburn) {

  load(sample_file)
  
  burn = 1:nburn
  
  # extract posterior samples for parameters to recover
  samples.names = colnames(samples)
  tgts = grep('(^pi\\[[13]\\]|^beta_lambda\\[[123], 1\\])', 
              samples.names, value = TRUE)
  m = mcmc(samples[-burn, tgts])
  
  # assemble posterior summaries
  df = data.frame(
    pt_est = colMeans(m),
    sd = apply(m, 2, sd),
    hpd = HPDinterval(m),
    truth = c((sim_pkg$consts$beta_lambda_prior_mean[, 'intercept']), 
              pi[c(1,3)]),
    tstep = as.numeric(str_extract(sample_file, '[0-9]+'))
  )
  
  df
  
}


