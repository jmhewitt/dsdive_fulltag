post_sim_targets = list(

  tar_target(
    name = covariate_combinations,
    command = expand.grid(
      prop_recent_deep = (1:12)/12,
      daytime = c(0,1),
      moonlit = c(0,1)
    )
  ),
  
  tar_target(
    name = deep_dive_time_preds_wrt_covariates, 
    command = {
      
      # location of MCMC files (also used as output path)
      # path = nim_fit
      path = 'output/mcmc/nim_fit/'
      
      # posterior parameter sample files
      param_sample_files = dir(
        path = path, pattern = 'parameter_samples_[0-9]+', full.names = TRUE
      )
      
      # column labels for posterior parameter samples
      param_label_file = dir(
        path = path, pattern = 'parameter_samples_column', full.names = TRUE
      )
      
      # load posterior parameter samples 
      param_samples = do.call(rbind, lapply(param_sample_files, function(f) {
        readRDS(f)
      }))
      
      # load data package
      nim_pkg = readRDS(
        file = dir(path = path, pattern = 'nim_pkg', full.names = TRUE)
      )
      
      # label posterior parameter samples
      colnames(param_samples) = readRDS(param_label_file)
      
      burn = 1:2e3
      
      # indices of posterior samples to process
      post_samples = (1:nrow(param_samples))[-burn]
      
      # subset posterior samples according to parallelization task
      post_samples = post_samples
      
      # build initial covariate matrix
      covariates = rbind(
        intercept = rep(1, 12),
        daytime = rep(covariate_combinations$daytime, 12), 
        moonlit = rep(covariate_combinations$moonlit, 12),
        prop_recent_deep = pmax(
          covariate_combinations$prop_recent_deep - (0:11)/12, 
          0
        ),
        prop_recent_deep3 = rep(0,12)
      )
      covariates['prop_recent_deep3',] = (
        covariates['prop_recent_deep',] - 0.5
      )^3
      
      # build initial depth vector
      depths = rep(x = c(16, 1), 
          c(covariate_combinations$prop_recent_deep * 12,
            12 - covariate_combinations$prop_recent_deep * 12)
      )
      
      # draw posterior predictive samples of time to next deep observation
      post_pred_samples = sapply(post_samples, function(sample_ind) {
        
        # extract stage transition coefficients
        betas_tx = array(dim = c(nim_pkg$consts$n_covariates,
                                 nim_pkg$consts$n_stages,
                                 nim_pkg$consts$n_stages - 1))
        for(i in 1:nrow(betas_tx)) {
          for(j in 1:ncol(betas_tx)) {
            for(k in 1:dim(betas_tx)[3]) {
              betas_tx[i,j,k] = param_samples[
                sample_ind, paste('beta_tx[', i, ', ', j, ', ', k, ']', 
                                  sep = '')
              ]
            }
          }
        }
        
        # forward-simulate dive from predictive distribution
        fwd_sim = fwd_sim_to_depth_fixed_covs(
          stages = rep(
            which(rownames(movement_classes$stage_defs) == 'slow_descent'), 
            12
          ),
          depths = depths, 
          covariates = covariates, 
          n_max = 1e3, 
          nim_pkg = nim_pkg,
          lambda = param_samples[sample_ind, c('lambda[1]', 'lambda[2]')], 
          betas_tx = betas_tx, 
          template_bins = template_bins,
          depth_threshold = deep_dive_depth,
          deeper = TRUE
        )
        
        # sampled time to deep depth  
        length(fwd_sim$depths) - 1
        
      })
      
      list(list(
        covariate_combination = covariate_combinations,
        samples = post_pred_samples
      ))
      
    }, 
    pattern = map(covariate_combinations), 
    deployment = 'worker',
    memory = 'transient'
  ),
  
  tar_target(
    name = deep_dive_time_interpretations_wrt_covariates,
    command = {
    
      deep_dive_time_preds_wrt_covariates = readRDS(
        'deep_dive_time_preds_wrt_covariates.rds'
      )
      
      # # long format of output, for plotting
      # df = do.call(
      #   rbind, lapply(deep_dive_time_preds_wrt_covariates, function(res) {
      #     
      #     tab = tabulate(res$samples)
      #     
      #     data.frame(
      #       res$covariate_combination,
      #       n = 1:length(tab),
      #       p = tab/sum(tab)
      #     ) %>% 
      #       dplyr::mutate(cdf = cumsum(p))
      # }))
      # 
      # # plot posterior predictive distributions
      # ggplot(df, aes(x = n * 5 / 60, y = p, col = factor(prop_recent_deep))) + 
      #   geom_line() + 
      #   facet_grid(daytime~moonlit) + 
      #   theme_few() + 
      #   xlim(0,3) + 
      #   ylim(0,.15)
      
      # summarize distributions (in hours units)
      df.summary = do.call(
        rbind, lapply(deep_dive_time_preds_wrt_covariates, function(res) {
          
          hpds = HPDinterval(mcmc(res$samples * 5 / 60))
          
          data.frame(
            res$covariate_combination,
            mean = mean(res$samples * 5 / 60),
            lwr = hpds[,'lower'],
            upr = hpds[,'upper']
          )
      }))
      
      # plot summaries of posterior predictive distributions
      group_sep = .015
      pl = ggplot(df.summary %>% 
               filter(
                 daytime * moonlit == 0,
                      prop_recent_deep < 1) %>% 
               mutate(
                 scenario = ifelse(daytime == 1, 'Daytime', 
                                   ifelse(moonlit == 1, 'Moonlit night', 
                                          'Dark night')),
                 scenario = factor(scenario, levels = c('Daytime', 
                                                        'Moonlit night', 
                                                        'Dark night')),
                 eps = ifelse(scenario == 'Daytime', group_sep, 
                              ifelse(scenario == 'Moonlit night', 0, 
                                     -group_sep))
               ), 
             aes(x = prop_recent_deep + eps, y = mean, ymin = lwr, ymax = upr,
                 col = scenario)) + 
        geom_pointrange(alpha = .8) + 
        scale_color_brewer('Celestial covariate', 
                           type = 'qual', palette = 'Dark2') + 
        scale_x_continuous(
          'Proportion of recent observations at depth', 
          breaks = c(0.08, .25, .5, .75, .92)
        ) + 
        ylab('Hours until next observed deep depth') + 
        theme_few() + 
        theme(panel.border = element_blank(), 
              panel.grid.major.y = element_line(colour = 'grey95'))
      
      # save plot
      f = file.path('output', 'covariate_effects')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      ggsave(pl, filename = file.path(f, paste(tar_name(), '.png', sep = '')), 
             dpi = 'print')
      
    })
  
)
