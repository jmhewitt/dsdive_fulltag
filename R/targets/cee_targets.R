cee_targets = list(
  
  tar_target(
    name = cee_predictive_samples,
    command = {
      
      # collate a posterior sample of parameters from across the multiple chains
      # that will be used to draw posterior predictive samples for covariate
      # effect interpretation, and cee response estimation
      
      #
      # enumerate random starts outputs
      #
      
      # scrape directory where posterior samples from chains live
      outdirs = dir(
        path = file.path('output', 'mcmc', 'fixed_init_beta'), 
        pattern = 'fit_marginalized_model_', 
        full.names = TRUE
      )
      
      # post-burn in samples to work with
      post_start = 8.5e3
      post_end = 10e3
      sample_inds = post_start:post_end
      
      # compute posterior summaries for each chain
      summaries = lapply(outdirs, function(d) {
        
        message(d)
        
        #
        # load, label, and merge posterior samples
        #
        
        mvSample_files = dir(
          path = d, 
          pattern = 'mvSamples_[0-9]+', 
          full.names = TRUE
        )
        
        if(length(mvSample_files) == 0) {
          return(NULL)
        }
        
        samples = do.call(rbind, lapply(mvSample_files, readRDS))
        
        # skip file if not enough samples are available for subsetting
        if(nrow(samples) < max(sample_inds)) {
          return(NULL)
        }
        
        colnames(samples) = readRDS(dir(
          path = d, 
          pattern = 'mvSamples_colnames',
          full.names = TRUE
        ))
        
        # target nodes
        tgt = c(
          paste('lambda[', 1:2, ']', sep = ''),
          paste('pi[', c(1,3), ']', sep = '')
        )
        
        # extract samples into single object
        list(
          chain = d,
          summaries = data.frame(
            param = tgt,
            mean = colMeans(samples[sample_inds, tgt]),
            sd = apply(samples[sample_inds, tgt], 2, sd),
            HPDinterval(mcmc(samples[sample_inds, tgt])),
            rep = d
          )
        )
      })
      
      # remove null output
      summaries = summaries[!sapply(summaries, is.null)]
      
      # # rules to identify chains that have converged to distant local modes
      # keep_chain = sapply(summaries, function(s) {
      #   if(s$summaries['pi[1]', 'mean'] > .9) {
      #     return(FALSE)
      #   }
      #   if(s$summaries['lambda[1]', 'mean'] > .7) {
      #     return(FALSE)
      #   }
      #   TRUE
      # })
      keep_chain = rep(TRUE, length(summaries))
      
      # look at convergence diagnostics
      df = do.call(
        rbind, lapply(summaries[keep_chain], function(s) s$summaries)
      ) %>% 
        mutate(rep = as.numeric(stringr::str_extract(rep, '[0-9]+')))
      pl = ggplot(df, aes(x = rep, y = mean, ymin = lower, ymax = upper)) + 
        geom_pointrange() + 
        facet_wrap(~param, labeller = label_parsed, scales = 'free_y') + 
        theme_few()
      print(pl)
      
      #
      # extract samples used for posterior predictive sampling
      #
      
      # chains to aggregate parameters from
      chains = sapply(summaries[keep_chain], function(s) s$chain)
      nchains = length(chains)
      
      # target number of aggregated samples to collect
      npost_samples = 1e4
      
      # number of samples to extract from each chain
      samples_per_chain = ceiling(npost_samples / nchains)
      
      # post-burn in posterior sample indices to extract from each chain
      sample_inds = seq(
        from = post_start, to = post_end, length.out = samples_per_chain
      )
      
      # aggregate parameter samples from each chain
      samples = do.call(rbind, lapply(chains, function(d) {
        
        message(d)
        
        #
        # load, label, and merge posterior samples
        #
        
        mvSample_files = dir(
          path = d, 
          pattern = 'mvSamples_[0-9]+', 
          full.names = TRUE
        )
        
        mvSample2_files = dir(
          path = d, 
          pattern = 'mvSamples2_[0-9]+', 
          full.names = TRUE
        )
        
        if(length(mvSample_files) == 0) {
          return(NULL)
        }
        
        if(length(mvSample2_files) == 0) {
          return(NULL)
        }
        
        samples = do.call(rbind, lapply(mvSample_files, readRDS))
        samples2 = do.call(rbind, lapply(mvSample2_files, readRDS))
        
        # skip file if not enough samples are available for subsetting
        if(nrow(samples) < max(sample_inds)) {
          return(NULL)
        }
        
        # skip file if not enough samples are available for subsetting
        if(nrow(samples2) < max(sample_inds)) {
          return(NULL)
        }
        
        colnames(samples) = readRDS(dir(
          path = d, 
          pattern = 'mvSamples_colnames',
          full.names = TRUE
        ))
        
        colnames(samples2) = readRDS(dir(
          path = d, 
          pattern = 'mvSamples2_colnames',
          full.names = TRUE
        ))
        
        # extract samples into single object
        cbind(samples, samples2)[sample_inds,]
      }))
      
      # save thinned, collated posterior samples 
      f = file.path('output', 'mcmc', 'fixed_init_beta')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      saveRDS(samples, file = file.path(f, paste(tar_name(), '.rds', sep = '')))
      
      0
    }
  ),

  tar_target(
    name = cee_results,
    command = {
      
      f = dir(path = file.path('output', 'cee', 'samples'), full.names = TRUE)
      
      cee_predictions = lapply(as.list(f), readRDS)
      
      cee_predictions_mod = lapply(cee_predictions, function(x) {
        if(!is.null(x[[1]]$pval)) {
          x[[1]]$observed_time_to_deep = c(
            x[[1]]$observed_time_to_deep,
            dive_interval = as.numeric(diff(x[[1]]$observed_time_to_deep))
          )
          samples = apply(x[[1]]$baseline_deep_pred_samples, 2, diff)
          x[[1]]$baseline_deep_pred_samples = rbind(
            x[[1]]$baseline_deep_pred_samples,
            interval = samples
          )
          x[[1]]$pval = c(
            x[[1]]$pval,
            interval_shorter = mean(
              samples < x[[1]]$observed_time_to_deep['dive_interval']
            ),
            interval_longer = mean(
              x[[1]]$observed_time_to_deep['dive_interval'] < samples
            )
          )
          x[[1]]$median = apply(x[[1]]$baseline_deep_pred_samples, 1, median)
          x[[1]]$sd = apply(x[[1]]$baseline_deep_pred_samples, 1, sd)
          x[[1]]$zscore = 
            (x[[1]]$observed_time_to_deep - 
             rowMeans(x[[1]]$baseline_deep_pred_samples)) / x[[1]]$sd
            
        }
        x
      })
      
      manifest = do.call(rbind, lapply(cee_predictions_mod, function(x) {
        cee_id = x[[1]]$cee
        if(inherits(cee_id, 'data.frame')) {
          cee_id = cee_id$cee_id
        }
        if(is.null(x[[1]]$pval)) {
          data.frame(
            tag = x[[1]]$tag,
            cee = cee_id,
            deep_at_exposure = NA,
            first_deep_shorter = NA,
            second_deep_shorter = NA,
            second_deep_longer = NA,
            interval_shorter = NA,
            interval_longer = NA,
            z_first_deep = NA,
            z_second_deep = NA,
            z_interval = NA,
            sd_first_deep = NA,
            sd_second_deep = NA,
            sd_interval = NA,
            median_first_deep = NA,
            median_second_deep = NA,
            median_interval = NA,
            obs_first_deep = NA,
            obs_second_deep = NA,
            obs_interval = NA,
            msg = x[[1]]$error
          )
        } else {
          data.frame(
            tag = x[[1]]$tag,
            cee = cee_id,
            deep_at_exposure = x[[1]]$deep_at_exposure,
            first_deep_shorter = x[[1]]$pval['first_deep_shorter'],
            second_deep_shorter = x[[1]]$pval['second_deep_shorter'],
            second_deep_longer = x[[1]]$pval['second_deep_longer'],
            interval_shorter = x[[1]]$pval['interval_shorter'],
            interval_longer = x[[1]]$pval['interval_longer'],
            z_first_deep = x[[1]]$zscore['first_deep'],
            z_second_deep = x[[1]]$zscore['second_deep'],
            z_interval = x[[1]]$zscore['dive_interval'],
            sd_first_deep = x[[1]]$sd['first_deep'],
            sd_second_deep = x[[1]]$sd['second_deep'],
            sd_interval = x[[1]]$sd['dive_interval'],
            median_first_deep = x[[1]]$median['first_deep'],
            median_second_deep = x[[1]]$median['second_deep'],
            median_interval = x[[1]]$median['interval'],
            obs_first_deep = x[[1]]$observed_time_to_deep['first_deep'],
            obs_second_deep = x[[1]]$observed_time_to_deep['second_deep'],
            obs_interval = x[[1]]$observed_time_to_deep['dive_interval'],
            msg = NA
          )
        }
      }))
      
      x = manifest %>% 
        filter(deep_at_exposure == FALSE) %>% 
        select(first_deep_shorter) %>% 
        unlist() %>% 
        as.numeric()
      
      x2f = manifest %>% 
        filter(deep_at_exposure == FALSE) %>% 
        select(interval_shorter) %>% 
        unlist() %>% 
        as.numeric()
      
      x2t = manifest %>% 
        filter(deep_at_exposure == TRUE) %>% 
        select(interval_shorter) %>% 
        unlist() %>% 
        as.numeric()
      
      # index 2 is also good!
      cee_results = cee_predictions_mod[[6]][[1]]
      e = ecdf(cee_results$baseline_deep_pred_samples['first_deep',]/3600)
      pl = ggplot() + 
        # posterior cdf
        stat_function(fun = function(x) e(x), n = 1e3) + 
        # posterior p-value
        geom_hline(
          yintercept = cee_results$pval['first_deep_shorter'],
          alpha = .1
        ) + 
        # observed value
        geom_vline(
          xintercept = cee_results$observed_time_to_deep['first_deep']/3600,
          lty = 3
        ) + 
        # axes
        scale_x_continuous(
          # 'Time to reach first deep depth (h)', 
          expression(H[1]('\u2113'[ij^'*'])),
          breaks = c(0,30,60,90,120)/60,
          limits = c(0, 2*60)/60
        ) + 
        scale_y_continuous(
          # 'Post. pred. CDF', 
          expression(F(H[1]('\u2113'[ij^'*']))),
          limits = c(0,1),
          sec.axis = sec_axis(
            trans = identity, 
            breaks = cee_results$pval['first_deep_shorter'], 
            # labels = paste(
            #   # 'F(x)=', 
            #   round(cee_results$pval['first_deep_shorter'],2),
            #   sep = ''
            # ),
            labels = expression(F(h[1]('\u2113'[ij^'*']))==.02)
          )
        ) +
        # title and formatting
        ggtitle(gsub(pattern = 'Tag', replacement = '', x = cee_results$tag)) + 
        theme_few()
      
      e2 = ecdf(cee_results$baseline_deep_pred_samples['interval',]/3600)
      pl2 = ggplot() + 
        # posterior cdf
        stat_function(fun = function(x) e2(x), n = 1e3) + 
        # posterior p-value
        geom_hline(
          yintercept = cee_results$pval['interval_shorter'],
          alpha = .1
        ) + 
        # observed value
        geom_vline(
          xintercept = cee_results$observed_time_to_deep['dive_interval']/3600,
          lty = 3
        ) + 
        # axes
        scale_x_continuous(
          # 'Deep dive recovery time (h)', 
          expression(H[2]('\u2113'[ij^'*'])-H[1]('\u2113'[ij^'*'])),
          breaks = seq(0,12,by = 2),
          limits = c(0,12)
        ) + 
        scale_y_continuous(
          # 'Post. pred. CDF',
          expression(F(H[2]('\u2113'[ij^'*'])-H[1]('\u2113'[ij^'*']))),
          limits = c(0,1),
          sec.axis = sec_axis(
            trans = identity, 
            breaks = cee_results$pval['interval_shorter'], 
            # labels = paste(
            #   'F(x)=', 
            #   round(cee_results$pval['interval_shorter'],2),
            #   sep = ''
            # )
            labels = expression(F(h[2]('\u2113'[ij^'*'])-h[1]('\u2113'[ij^'*']))==1)
          )
        ) +
        # title and formatting
        ggtitle(gsub(pattern = 'Tag', replacement = '', x = cee_results$tag)) + 
        theme_few()
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      ggsave(
        ggarrange(
          pl + theme(plot.title = element_blank()), 
          pl2 + theme(plot.title = element_blank()), 
                      # axis.title.y = element_blank()),
          nrow = 1, ncol = 2, labels = 'AUTO'
          # align = 'v'
          # label.x = c(0,-.05)
        ),
        filename = file.path(f, paste(tar_name(), '_', cee_results$tag, '.png',
                                      sep = '')),
        width = 8, height = 3, dpi = 'print'
      )
      
      # qq-uniform plot; the model probs. are on the x-axis, and their sample
      # percentile is on the y-axis.  If the animals were not impacted, then 
      # the percentiles would be "balanced", and lie closer to the 1:1 line.
      # # We don't see this, in fact, we see a shift toward smaller probs.
      # plo = ggarrange(
      #   
      #   ggplot(data.frame(x), aes(x=x)) + 
      #     geom_vline(xintercept = 0.5, col = 'grey90') +
      #     stat_ecdf(geom = 'point') + 
      #     xlab(expression(t)) + 
      #     ylab(expression('Freq. '~bgroup('{',F(bold('\u00B7'))<t,'}'))) +
      #     # ylab(expression('Freq. '~bgroup('{',F(symbol('\267'))<t,'}'))) + 
      #     geom_abline(slope = 1, intercept = 0, lty = 3) + 
      #     # ggtitle('First deep depth percentiles') + 
      #     ggtitle(expression(F(H[1]('\u2113'[ij^'*']))~' comparison'))+
      #     theme_few() + 
      #     theme(
      #       plot.title = element_text(hjust = .5),
      #       axis.title.y = element_text(angle = 0, vjust = .5)) + 
      #     xlim(0,1),
      #   
      #   ggplot(data.frame(x2f), aes(x=x2f)) + 
      #     geom_vline(xintercept = 0.5, col = 'grey90') +
      #     stat_ecdf(geom = 'point') + 
      #     xlab(expression(t)) + 
      #     ylab(expression('Freq. '~bgroup('{',F(bold('\u00B7'))<t,'}'))) +
      #     # ylab(expression('Freq. '~bgroup('{',F(symbol('\267'))<t,'}'))) + 
      #     geom_abline(slope = 1, intercept = 0, lty = 3) + 
      #     # ggtitle('Deep dive recovery percentiles') + 
      #     ggtitle(expression(
      #       F(H[2]('\u2113'[ij^'*'])-H[1]('\u2113'[ij^'*']))~' comparison'
      #     ))+
      #     theme_few() +
      #     theme(plot.title = element_text(hjust = .5),
      #           axis.title.y = element_blank()) + 
      #     xlim(0,1),
      #   
      #   nrow = 1, ncol = 2, labels = 'AUTO', align = 'v'
      # )
      
      # facet version
      # repeat axis title solution: https://stackoverflow.com/questions/70022117/how-to-add-y-axis-title-for-each-facet-row-in-ggplot
      plo = rbind(
        data.frame(x = x, series = "F(H[1]('\u2113'[ij^'*']))~' comparison'"),
        data.frame(x = x2f, series = "F(H[2]('\u2113'[ij^'*'])-H[1]('\u2113'[ij^'*']))~' comparison'")
      ) %>% ggplot(aes(x=x)) + 
        geom_vline(xintercept = .5, col = 'grey90') + 
        stat_ecdf(geom = 'point') +
        geom_abline(slope = 1, intercept = 0, lty = 3) + 
        xlab(expression(t)) + 
        ylab(expression('Freq. '~bgroup('{',F(bold('\u00B7'))<t,'}'))) +
        facet_wrap(~series, labeller = label_parsed, scales = 'free_x') +
        theme_few() + 
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
        xlim(0,1)
      
      gt = ggplot_gtable(ggplot_build(plo))
      which.xlab = grep('xlab-b', gt$layout$name)
      gt = gtable::gtable_add_grob(gt, gt$grobs[which.xlab], 10, 5)
      gt = gtable::gtable_add_grob(gt, gt$grobs[which.xlab], 10, 9)
      # remove the original axis title
      gt = gtable::gtable_filter(gt, 'xlab-b', invert = TRUE) 

      ggsave(gt, 
             filename = file.path(f, paste(tar_name(), '_pop.png', sep ='')),
             width = 8, height = 3, dpi = 'print')
      
      ind = 10
      ggplot(
        data.frame(t(cee_predictions[[ind]][[1]]$baseline_deep_pred_samples)),
        aes( x = first_deep/3600, y = (second_deep-first_deep)/3600)
      ) + 
        geom_point(alpha = .05) + 
        geom_point(
          data = data.frame(t(cee_predictions[[ind]][[1]]$observed_time_to_deep)),
          col = 'green'
        ) + 
        stat_density_2d() +
        ggtitle(cee_predictions[[ind]][[1]]$tag) + 
        theme_few()
      
      xlsx::write.xlsx(
        manifest, 'DiveProbs.xlsx', row.names = FALSE, showNA = FALSE
      )
      
    }
  )  
)
