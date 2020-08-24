priors_and_inits = function(nim_data, template_bins, sattag_timestep, 
                            expm_delta) {

  # initialize output
  nim_pkg = nim_data
  nim_pkg$inits = list()
  
  # add timestep and approximation, for likelihood computation
  nim_pkg$consts$tstep = sattag_timestep
  nim_pkg$consts$delta = 1e-10

  #
  # priors for dive durations
  #

  # reformat dive data for 85% rule function
  dives.obs = apply(nim_pkg$consts$dive_relations, 1, function(r) {
    list(
      depth.bins = template_bins,
      dive = list(
        depths = nim_pkg$data$depths[r['depth_first']:r['depth_last']],
        times = nim_pkg$data$times[r['depth_first']:r['depth_last']] -
          nim_pkg$data$times[r['depth_first']]
      )
    )
  })

  # estimate stage durations by tag
  durations.est = cbind(
    times.stages(dives.obs) * 60,
    tag = nim_pkg$consts$dive_relations[,'tag']
  )

  colnames(durations.est) = c('sub.time.sec', 'bottom.time.sec', 'tag')

  # stage duration priors based on log-normal marginals
  duration_priors = do.call(
    rbind, lapply(sort(unique(durations.est$tag)), function(i) {
      # identify estimates for ith tag
      inds = durations.est$tag == i
      # estimate dive durations and correlation on log scale
      res = c(
        fitdistr(durations.est$sub.time.sec[inds], 'lognormal')$estimate,
        fitdistr(durations.est$bottom.time.sec[inds], 'lognormal')$estimate,
        cov(log(durations.est$sub.time.sec[inds]),
            log(durations.est$bottom.time.sec[inds]))
      )
      # format and return
      names(res) = c('G1_mean', 'G1_sd', 'G2_mean', 'G2_sd', 'cov_log')
      res
    }))
  
  # stage duration priors based on gamma marginals
  duration_priors_gamma = do.call(
    rbind, lapply(sort(unique(durations.est$tag)), function(i) {
      # identify estimates for ith tag
      inds = durations.est$tag == i
      # estimate dive durations and correlation on log scale
      res = c(
        fitdistr(durations.est$sub.time.sec[inds], 'gamma')$estimate,
        fitdistr(durations.est$bottom.time.sec[inds], 'gamma')$estimate,
        cor(log(durations.est$sub.time.sec[inds]),
            log(durations.est$bottom.time.sec[inds]))
      )
      # format and return
      names(res) = c('G1_shape', 'G1_rate', 'G2_shape', 'G2_rate', 'cor_log')
      res
    }))
  
  # prior means for stage durations
  lognormal_means = apply(duration_priors, 1, function(r) {
    c(exp(r['G1_mean'] + r['G1_sd']^2/2),
      exp(r['G2_mean'] + r['G2_sd']^2/2))
  })
  
  # prior means for stage durations
  gamma_means = apply(duration_priors_gamma, 1, function(r) {
    c(r['G1_shape']/r['G1_rate'],
      r['G2_shape']/r['G2_rate'])
  })
  
  # prior sd's for stage durations
  lognormal_sd = apply(duration_priors, 1, function(r) {
    sqrt(c( (exp(r['G1_sd']^2) - 1) * exp(2 * r['G1_mean'] + r['G1_sd']^2),
            (exp(r['G2_sd']^2) - 1) * exp(2 * r['G2_mean'] + r['G2_sd']^2)
      ))
  })
  
  # prior sd's for stage durations
  gamma_sd = apply(duration_priors_gamma, 1, function(r) {
    sqrt(c(r['G1_shape']/r['G1_rate']^2,
           r['G2_shape']/r['G2_rate']^2))
  })
  
  # lognormal vs gamma means only differ by a few percent
  round((lognormal_means[2,]/gamma_means[2,] - 1) * 100)
  
  # lognormal sd's tend to be greater by order of magnitude in pct
  round((lognormal_sd[2,]/gamma_sd[2,] - 1) * 100)
  
  # this yields an average of 1.5 minutes wrt sd
  summary((lognormal_sd[2,] - gamma_sd[2,])/60)
  
  # lognormal parameters that yield equivalent means and sd's to the gamma fits
  alt_params = cbind(
    matrix(
      lognormal.params(mean = gamma_means, var = gamma_sd^2), ncol = 2
    )[1:nim_pkg$consts$N_tags, , drop = FALSE],
    matrix(
      lognormal.params(mean = gamma_means, var = gamma_sd^2), ncol = 2
    )[-(1:nim_pkg$consts$N_tags), , drop = FALSE],
    NA
  )
  alt_params[,c(2,4)] = sqrt(alt_params[,c(2,4)])
  alt_params[, 5] = duration_priors_gamma[, 'cor_log'] * 
    alt_params[, 2] * alt_params[, 4]
  colnames(alt_params) = colnames(duration_priors)
  
  # # effect of alternate priors is largely to tighten up the distribution tails
  # lapply(7:1, function(i) {
  #   inds = durations.est$tag == i
  #   
  #   par(mfrow = c(1,1))
  #   
  #   plot(density((durations.est$bottom.time.sec[inds])), main = i)
  #   curve(dlnorm(x, 
  #               meanlog = alt_params[i, 'G2_mean'], 
  #               sdlog = alt_params[i, 'G2_sd']), col = 2, add = TRUE)
  #   curve(dlnorm(x, 
  #               meanlog = duration_priors[i, 'G2_mean'], 
  #               sdlog = duration_priors[i, 'G2_sd']), col = 4, add = TRUE)
  #   curve(dgamma(x, 
  #                shape = duration_priors_gamma[i, 'G2_shape'],
  #                rate = duration_priors_gamma[i, 'G2_rate']), col = 3, add = TRUE)
  #   
  #   par(mfrow = c(1,1))
  #   
  # })
  
  # use alternate priors
  duration_priors = alt_params

  # reformat duration priors for nimble
  nim_pkg$consts$xi_prior_means = duration_priors[,c('G1_mean', 'G2_mean'), 
                                                  drop = FALSE]
  nim_pkg$consts$xi_prior_covs = array(data = NA,
                                       dim = c(nrow(duration_priors), 2, 2))
  for(i in 1:nrow(duration_priors)) {
    nim_pkg$consts$xi_prior_covs[i,1,1] = duration_priors[i, 'G1_sd']^2
    nim_pkg$consts$xi_prior_covs[i,1,2] = duration_priors[i, 'cov_log']
    nim_pkg$consts$xi_prior_covs[i,2,1] = duration_priors[i, 'cov_log']
    nim_pkg$consts$xi_prior_covs[i,2,2] = duration_priors[i, 'G2_sd']^2
  }


  #
  # random effect and covariate priors
  #

  # prior means for logit-pi model coefficients
  nim_pkg$consts$pi_prior = rbind(
    c(shape1 = 5, shape2 = 2),
    c(shape1 = 1, shape2 = 1),
    c(shape1 = 2, shape2 = 5)
  )

  # prior means for log-lambda model coefficient
  nim_pkg$consts$beta_lambda_prior_mean = rbind(
    c(intercept = log(1.5), sex = 0),
    c(intercept = log(.3), sex = 0),
    c(intercept = log(1), sex = 0)
  )

  # prior uncertainty for log-lambda model coefficients
  nim_pkg$consts$beta_lambda_prior_sd = rbind(
    c(intercept = 1e2, sex = 1e2),
    c(intercept = 1e2, sex = 1e2),
    c(intercept = 1e2, sex = 1e2)
  )

  # prior distribution for log-lambda model random effect scale
  nim_pkg$consts$sigma_lambda_priors = rbind(
    c(shape = 2, rate = 1),
    c(shape = 2, rate = 1),
    c(shape = 2, rate = 1)
  )


  #
  # set initial values for random variables
  #

  nim_pkg$inits$endpoints = rowMeans(nim_pkg$consts$endpoint_priors)

  nim_pkg$inits$xi = as.matrix(
    durations.est[,c('sub.time.sec', 'bottom.time.sec')]
  )
  colnames(nim_pkg$inits$xi) = c('sub_time_sec', 'bottom_time_sec')

  nim_pkg$inits$log_xi = log(nim_pkg$inits$xi)

  nim_pkg$inits$T = cbind(
    T0 = nim_pkg$inits$endpoints[nim_pkg$consts$dive_relations[,'T0_endpoint']],
    T1 = NA,
    T2 = NA,
    T3 = nim_pkg$inits$endpoints[nim_pkg$consts$dive_relations[,'T3_endpoint']]
  )

  nim_pkg$inits$T[,'T1'] = nim_pkg$inits$T[,'T0'] +
    nim_pkg$inits$xi[,'sub_time_sec']

  nim_pkg$inits$T[,'T2'] = nim_pkg$inits$T[,'T1'] +
    nim_pkg$inits$xi[,'bottom_time_sec']

  # verify all time random effects are well ordered
  stopifnot(all(
    nim_pkg$inits$T[,'T0'] < nim_pkg$inits$T[,'T1'],
    nim_pkg$inits$T[,'T1'] < nim_pkg$inits$T[,'T2'],
    nim_pkg$inits$T[,'T2'] < nim_pkg$inits$T[,'T3']
  ))
  
  nim_pkg$inits$pi = c(.99, .5, .01)
  
  nim_pkg$inits$logit_pi = qlogis(nim_pkg$inits$pi)

  nim_pkg$inits$beta_lambda = nim_pkg$consts$beta_lambda_prior_mean

  nim_pkg$inits$log_lambda = matrix(
    rep(nim_pkg$consts$beta_lambda_prior_mean[,'intercept'],
        nim_pkg$consts$N_tags),
    ncol = 3, byrow = TRUE
  )

  nim_pkg$inits$lambda = exp(nim_pkg$inits$log_lambda)

  nim_pkg$inits$sigma_lambda = rep(1, 3)
  
  saveRDS(nim_pkg, paste(id_chr(), '.rds', sep = ''))
  
  nim_pkg

}
