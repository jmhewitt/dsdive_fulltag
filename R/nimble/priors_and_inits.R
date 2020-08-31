priors_and_inits = function(nim_data, template_bins, sattag_timestep, 
                            expm_delta, mcmc_sample_dir) {

  # initialize output
  nim_pkg = readRDS(nim_data)
  nim_pkg$inits = list()
  
  # add timestep and approximation, for likelihood computation
  nim_pkg$consts$tstep = sattag_timestep
  nim_pkg$consts$delta = 1e-10
  
  
  #
  # priors and inits for stage duration parameters
  #
  
  # deep dives
  nim_pkg$consts$G_prior_mean = c(0, 0)
  nim_pkg$consts$G_prior_sd = c(1e2, 1e2)
  nim_pkg$consts$G_prior_shape = c(2, 2)
  nim_pkg$consts$G_prior_rate = c(1, 1)
  nim_pkg$consts$G_prior_cor = c(1, 1)
  # nim_pkg$inits$xi_prior_cor_scaled = .25
  nim_pkg$inits$xi_prior_cor_scaled = .8
  
  # shallow dives
  nim_pkg$consts$G_shallow_prior_mean = c(0, 1e2)
  nim_pkg$consts$G_shallow_prior_sd = c(2, 1)
  
  
  #
  # priors for dive durations
  #

  # reformat deep dive data for 85% rule function
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

  # estimate deep dive stage durations by tag
  durations.est = cbind(
    times.stages(dives.obs) * 60,
    tag = nim_pkg$consts$dive_relations[,'tag']
  )

  colnames(durations.est) = c('sub.time.sec', 'bottom.time.sec', 'tag')
  
  # estimate shallow dive descent duration as half of dive
  durations.shallow.est = data.frame(
    apply(
      nim_pkg$consts$dive_relations_shallow, 1, function(r) {
        diff(c(
          mean(nim_pkg$consts$endpoint_priors[r['T0_endpoint'],]),
          mean(nim_pkg$consts$endpoint_priors[r['T2_endpoint'],])
        ))/2
    }),
    tag = nim_pkg$consts$dive_relations_shallow[,'tag']
  )
  
  colnames(durations.shallow.est)[1] = 'sub.time.sec'
  
  # shallow dive stage duration priors based on gamma marginals
  duration_priors_shallow = do.call(
    rbind, lapply(sort(unique(durations.est$tag)), function(i) {
      # identify estimates for ith tag
      inds = durations.shallow.est$tag == i
      # develop empirical prior
      res = lognormal_moment_fit(x = durations.shallow.est$sub.time.sec[inds])
      # format and return
      names(res) = c('G1_mean', 'G1_var')
      res
    })
  )
  
  # deep dive stage duration priors based on gamma marginals
  duration_priors = do.call(
    rbind, lapply(sort(unique(durations.est$tag)), function(i) {
      # identify estimates for ith tag
      inds = durations.est$tag == i
      # develop empirical prior
      res = c(
        lognormal_moment_fit(x = durations.est$sub.time.sec[inds]),
        lognormal_moment_fit(x = durations.est$bottom.time.sec[inds]),
        cor(log(durations.est$sub.time.sec[inds]),
            log(durations.est$bottom.time.sec[inds]))
      )
      # format and return
      names(res) = c('G1_mean', 'G1_var', 'G2_mean', 'G2_var', 'cor_log')
      res
    })
  )
  
  # convert correlation to covariance
  duration_priors[, 'cor_log'] = duration_priors[, 'cor_log'] * 
    sqrt(duration_priors[, 'G1_var']) * sqrt(duration_priors[, 'G2_var'])
  colnames(duration_priors)[5] = 'cov_log'
  
  # reformat deep dive duration priors for nimble
  # nim_pkg$inits$xi_prior_means = duration_priors[,c('G1_mean', 'G2_mean'), 
  #                                                 drop = FALSE]
  nim_pkg$inits$xi_prior_means = matrix(0, nrow = 1, ncol = 2)
  
  nim_pkg$inits$xi_prior_covs = array(data = 1,
                                       dim = c(nrow(duration_priors), 2, 2))
  
  # for(i in 1:nrow(duration_priors)) {
  #   nim_pkg$inits$xi_prior_covs[i,1,1] = duration_priors[i, 'G1_var']
  #   nim_pkg$inits$xi_prior_covs[i,1,2] = duration_priors[i, 'cov_log']
  #   nim_pkg$inits$xi_prior_covs[i,2,1] = duration_priors[i, 'cov_log']
  #   nim_pkg$inits$xi_prior_covs[i,2,2] = duration_priors[i, 'G2_var']
  # }
  # 
  # reformat shallow dive duration priors for nimble
  # nim_pkg$inits$xi_shallow_prior = duration_priors_shallow
  nim_pkg$inits$xi_shallow_prior = matrix(0, nrow = 1, ncol = 2)


  #
  # random effect and covariate priors
  #
  
  # prior means for logit-pi model coefficients
  nim_pkg$consts$pi_prior = rbind(
    # deep dives
    c(shape1 = 5, shape2 = 2),
    c(shape1 = 1, shape2 = 1),
    c(shape1 = 2, shape2 = 5),
    # shallow dives
    c(shape1 = 5, shape2 = 2),
    c(shape1 = 2, shape2 = 5)
  )

  # prior means for log-lambda model coefficient
  nim_pkg$consts$beta_lambda_prior_mean = rbind(
    # deep dives
    c(intercept = log(1.5), sex = 0),
    c(intercept = log(.3), sex = 0),
    c(intercept = log(1), sex = 0),
    # shallow dives
    c(intercept = log(.6), sex = 0),
    c(intercept = log(.6), sex = 0)
  )

  # prior uncertainty for log-lambda model coefficients
  nim_pkg$consts$beta_lambda_prior_sd = rbind(
    # deep dives
    c(intercept = 1e2, sex = 1e2),
    c(intercept = 1e2, sex = 1e2),
    c(intercept = 1e2, sex = 1e2),
    # shallow dives
    c(intercept = 1e2, sex = 1e2),
    c(intercept = 1e2, sex = 1e2)
  )

  # prior distribution for log-lambda model random effect scale
  nim_pkg$consts$sigma_lambda_priors = rbind(
    # deep dives
    c(shape = 2, rate = 1),
    c(shape = 2, rate = 1),
    c(shape = 2, rate = 1),
    # shallow dives
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
  
  nim_pkg$inits$xi_shallow = as.numeric(durations.shallow.est$sub.time.sec)

  nim_pkg$inits$log_xi = log(nim_pkg$inits$xi)
  
  nim_pkg$inits$log_xi_shallow = log(nim_pkg$inits$xi_shallow)

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
  
  nim_pkg$inits$T_shallow = cbind(
    T0 = nim_pkg$inits$endpoints[
      nim_pkg$consts$dive_relations_shallow[,'T0_endpoint']],
    T1 = NA,
    T2 = nim_pkg$inits$endpoints[
      nim_pkg$consts$dive_relations_shallow[,'T2_endpoint']]
  )
  
  nim_pkg$inits$T_shallow[,'T1'] = nim_pkg$inits$T_shallow[,'T0'] +
    nim_pkg$inits$xi_shallow
  
  nim_pkg$inits$pi = c(.99, .5, .01, .99, .01)
  
  nim_pkg$inits$logit_pi = qlogis(nim_pkg$inits$pi)
  
  nim_pkg$inits$beta_lambda = nim_pkg$consts$beta_lambda_prior_mean

  nim_pkg$inits$log_lambda = matrix(
    rep(nim_pkg$consts$beta_lambda_prior_mean[,'intercept'],
        nim_pkg$consts$N_tags),
    ncol = 5, byrow = TRUE
  )

  nim_pkg$inits$lambda = exp(nim_pkg$inits$log_lambda)

  nim_pkg$inits$sigma_lambda = rep(1, 5)
  
  # save package
  dir.create(mcmc_sample_dir, recursive = TRUE, showWarnings = FALSE)
  f = file.path(mcmc_sample_dir, paste(id_chr(), '.rds', sep = ''))
  saveRDS(nim_pkg, f)
  
  # return path
  f
}
