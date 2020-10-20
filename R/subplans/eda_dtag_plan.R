eda_dtag_plan = drake_plan(

  # output directory for eda files
  eda_dtag_dir = {
    p = file.path('output', 'eda', 'dtag')
    dir.create(path = p, recursive = TRUE, showWarnings = FALSE)
    p
  },
  
  # output directory for eda plots
  eda_dtag_plots_dir = {
    p = file.path('plots', 'eda', 'dtag')
    dir.create(path = p, recursive = TRUE, showWarnings = FALSE)
    p
  },
  
  # data, in list format
  standardized_dtags = target(
    command = {
      dat = load_raw_dtag(dtag_files = dtag_files, cee_starts = cee_starts)
      f = file.path(eda_dtag_dir, paste(id_chr(), '.rds', sep =''))
      saveRDS(dat[[1]], file = f)
      f
    },
    dynamic = map(dtag_files),
    format = 'file'
  ),
  
  # Monte Carlo pre/post sample pairs
  pre_post_pairs_dtag = target(
    command = {
      # generate samples
      res = lapply(standardized_dtags, function(tag) {
        # load tag
        f = tag
        tag = readRDS(f)
        # package samples
        list(
          # tag ID
          tag = tag$tag,
          # configuration
          response_lag = response_lag,
          window_length = window_length,
          # samples
          pairs = sample_pairs(data = function() { data.frame(readRDS(f)) },
                               times = tag$times, 
                               exposure_time = tag$exposure_time,
                               response_lag = response_lag, 
                               window_length = window_length,
                               nsamples = 1e3, 
                               sparse = TRUE)
        )
      })
      # save samples, and return file name
      combo = paste(res[[1]]$tag, 
                    paste('lag', response_lag, sep = ''),
                    paste('length', window_length, sep = ''),
                    sep = '_')
      f = file.path(eda_dtag_dir, paste(id_chr(), '_', combo, '.rds', 
                                   sep =''))
      saveRDS(res, file = f)
      f
    },
    dynamic = map(standardized_dtags), 
    transform = cross(response_lag = c(0, .5, 1),
                      window_length = c(.125, .25, .5, 1, 2, 3, 5)),
    format = 'file'
  ),
  
  # evaluate summary statistics and MC p-value for pre/post pairs
  pre_post_tests_dtag = target(
    command = {
      # get function from string
      test_fn = getFunction(test_fn_name)
      # load sampled and observed pre/post pairs
      mcout = readRDS(pre_post_pairs_dtag)[[1]]
      # initialize return
      res = mcout
      res$pairs = NULL
      res$stat = test_fn_name
      # load data if needed
      if(isTRUE(mcout$pairs$observed$sparse)) {
        dat = mcout$pairs$observed$data_loader()
        mcout$pairs$observed$pre = dat[
          mcout$pairs$observed$pre[1]:mcout$pairs$observed$pre[2], , 
          drop = FALSE
        ]
        # expand sample
        mcout$pairs$observed$post = dat[
          mcout$pairs$observed$post[1]:mcout$pairs$observed$post[2], , 
          drop = FALSE
          ]
      }
      if(!is.null(mcout$pairs)) {
        if(length(mcout$pairs$resampled) > 0) {
          # evaluate test statistic on all pre/post samples
          res$resampled = sapply(mcout$pairs$resampled, function(sample) {
            if(is.null(sample)) {
              NA
            } else {
              # expand sample if needed
              if(isTRUE(sample$sparse)) {
                sample$pre = dat[sample$pre[1]:sample$pre[2], , drop = FALSE]
                sample$post = dat[sample$post[1]:sample$post[2], , drop = FALSE]
              }
              # compute summary statistic
              test_fn(sample)
            }
          })
          # compute test statistic and Monte Carlo p-value
          res$observed = test_fn(mcout$pairs$observed)
          res$pval = mean(res$resampled >= res$observed, na.rm = TRUE)
        }
      }
      # save output, and return file name
      combo = paste(res$tag,
                    paste('lag', res[['response_lag']], sep = ''),
                    paste('length', res[['window_length']], sep = ''),
                    sep = '_')
      f = file.path(eda_dir, paste(id_chr(), '_', combo, '.rds', sep =''))
      saveRDS(res, file = f)
      f
    }, 
    dynamic = map(pre_post_pairs_dtag),
    transform = cross(
      pre_post_pairs_dtag, 
      test_fn_name = c('total_absolute_dtag', 'mean_upward_velocity',
                       'mean_downward_velocity')
    ),
    format = 'file'
  ),
  
  # aggregate/organize lists of pre_post files by tagId
  pre_post_files_dtag = target(
    command = list(dplyr::combine(pre_post_tests_dtag)),
    dynamic = map(pre_post_tests_dtag),
    transform = combine(pre_post_tests_dtag)
  ),
  
  eda_tagsummary_dtag = target(
    command = {
      
      tag_files = pre_post_files_dtag[[1]]
      
      # name of the tag being processed
      targetTag = readRDS(tag_files[1])$tag
      
      # extract p-values for MC tests
      df = do.call(rbind, lapply(tag_files, function(f) {
        mcout = readRDS(f)
        data.frame(response_lag = mcout$response_lag,
                   window_length = mcout$window_length,
                   stat = mcout$stat,
                   pval = ifelse(is.null(mcout$pval), NA, mcout$pval))
      }))
      
      # plot output
      pl = ggplot(data = df %>%
                    mutate(signif = cut(as.numeric(pval),
                                        breaks = c(0, .05, .1, 1)),
                           window_length = factor(window_length),
                           response_lag = factor(
                             x = paste(response_lag, 'h lag', sep = ''),
                             levels = paste(sort(unique(response_lag)),
                                            'h lag', sep ='')
                           )
                    ),
                  mapping = aes(x = window_length, y = stat, fill = signif)) +
        facet_wrap(~response_lag) +
        geom_tile(height = .9, width = .9, col = 1, alpha = .6) +
        scale_fill_brewer('p val.', type = 'seq', palette = 'PuRd',
                          direction = -1) +
        ylab('Statistic') +
        xlab('Window length (h)') +
        coord_equal(ratio = .75) +
        theme_few() +
        ggtitle(targetTag)
      
      # save output
      f = file.path(eda_dtag_plots_dir,
                    paste(targetTag, '_', id_chr(), '.png', sep = ''))
      ggsave(pl, filename = f)
      f
    },
    dynamic = map(pre_post_files_dtag),
    format = 'file'
  ),
  
  dtag_eda_info = target(
    command = {
      
    }
  ),
  
  dtag_eda_plots = target(
    command = {
      # amount of time (min) to plot on either side of windows
      buffer = 0.5 * 60
      
      # load sampled and observed pre/post pairs
      mcout = readRDS(pre_post_pairs_dtag)[[1]]
      
      # load data if needed
      if(isTRUE(mcout$pairs$observed$sparse)) {
        dat = mcout$pairs$observed$data_loader()
      }
      
      # temporal extents of pre/post windows
      analysis_windows = data.frame(
        xmin = c(mcout$pairs$windows$pre[1], mcout$pairs$windows$post[1]),
        xmax = c(mcout$pairs$windows$pre[2], mcout$pairs$windows$post[2])
      )
      
      # plot data
      pl = ggplot(
        data = dat %>% 
          filter(mcout$pairs$windows$pre[1] - minutes(buffer) <= times, 
                 times <= mcout$pairs$windows$post[2] + minutes(buffer)), 
        aes(x = times, y = p)
      ) + 
        # base data
        geom_line() + 
        # pre/post-exposure windows
        geom_rect(data = analysis_windows,
                  mapping = aes(xmin = xmin, xmax = xmax,
                                ymin = 0, ymax = Inf), inherit.aes = FALSE,
                  fill = NA, col = 'black', lwd = 2) +
        # exposure time
        geom_vline(xintercept = dat$exposure_time[1], lty = 1, col = 'yellow') + 
        # scales
        scale_y_reverse('Depth (m)') + 
        scale_x_datetime('Time') + 
        # formatting
        theme_few() + 
        theme(panel.border = element_blank())
      
      # save output, and return file name
      combo = paste(mcout$tag,
                    paste('lag', mcout[['response_lag']], sep = ''),
                    paste('length', mcout[['window_length']], sep = ''),
                    sep = '_')
      f = file.path(eda_dtag_plots_dir, 
                    paste(id_chr(), '_', combo, sep =''))
      for(extension in c('.pdf', '.png')) {
        ggsave(pl, filename = paste(f, extension, sep =''))
      }
      f
    }, 
    transform = map(pre_post_pairs_dtag),
    format = 'file'
  )
  

)