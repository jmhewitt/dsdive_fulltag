cee_targets = list(

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
          'Time to reach first deep depth (h)', 
          breaks = c(0,30,60,90,120)/60,
          limits = c(0, 2*60)/60
        ) + 
        scale_y_continuous(
          'Post. pred. CDF', 
          limits = c(0,1),
          sec.axis = sec_axis(
            trans = identity, 
            breaks = cee_results$pval['first_deep_shorter'], 
            labels = paste(
              'F(x)=', 
              round(cee_results$pval['first_deep_shorter'],2),
              sep = ''
            )
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
          'Deep dive recovery time (h)', 
          breaks = seq(0,8,by = 2),
          limits = c(0,8)
        ) + 
        scale_y_continuous(
          'Post. pred. CDF',
          limits = c(0,1),
          sec.axis = sec_axis(
            trans = identity, 
            breaks = cee_results$pval['interval_shorter'], 
            labels = paste(
              'F(x)=', 
              round(cee_results$pval['interval_shorter'],2),
              sep = ''
            )
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
          pl2 + theme(plot.title = element_blank(), 
                      axis.title.y = element_blank()),
          nrow = 1, ncol =2, labels = 'AUTO', 
          label.x = c(0,-.05)
        ),
        filename = file.path(f, paste(tar_name(), '_', cee_results$tag, '.pdf',
                                      sep = '')),
        width = 8, height = 3
      )
      
      
      # qq-uniform plot; the model probs. are on the x-axis, and their sample
      # percentile is on the y-axis.  If the animals were not impacted, then 
      # the percentiles would be "balanced", and lie closer to the 1:1 line.
      # We don't see this, in fact, we see a shift toward smaller probs.
      plo = ggarrange(
        
        ggplot(data.frame(x), aes(x=x)) + 
          geom_vline(xintercept = 0.5, col = 'grey90') +
          stat_ecdf(geom = 'point') + 
          xlab('Post. predictive percentile') + 
          ylab('Empirical percentile') + 
          geom_abline(slope = 1, intercept = 0, lty = 3) + 
          ggtitle('First deep depth percentiles')+ 
          theme_few() + 
          theme(plot.title = element_text(hjust = .5)) + 
          xlim(0,1),
        
        ggplot(data.frame(x2f), aes(x=x2f)) + 
          geom_vline(xintercept = 0.5, col = 'grey90') +
          stat_ecdf(geom = 'point') + 
          xlab('Post. predictive percentile') + 
          ylab('Empirical percentile') + 
          geom_abline(slope = 1, intercept = 0, lty = 3) + 
          ggtitle('Deep dive recovery percentiles') + 
          theme_few() +
          theme(plot.title = element_text(hjust = .5),
                axis.title.y = element_blank()) + 
          xlim(0,1),
        
        nrow = 1, ncol = 2, labels = 'AUTO'
      )
      
      ggsave(plo, 
             filename = file.path(f, paste(tar_name(), '_pop.pdf', sep ='')),
             width = 8, height = 3)
      
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
