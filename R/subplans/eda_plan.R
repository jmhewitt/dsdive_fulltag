eda_plan = drake_plan(
  
  # output directory for eda files
  eda_dir = {
    p = file.path('output', 'eda')
    dir.create(path = p, recursive = TRUE, showWarnings = FALSE)
    p
  },
  
  # output directory for eda plots
  eda_plots_dir = {
    p = file.path('plots', 'eda')
    dir.create(path = p, recursive = TRUE, showWarnings = FALSE)
    p
  },
  
  # data, in list format
  standardized_sattags = target(
    load_raw(depth_files = depth_files, template_bins = template_bins, 
             cee_starts = cee_starts, dive_labels = dive_labels),
    dynamic = map(depth_files)
  ),
  
  # Monte Carlo pre/post sample pairs
  pre_post_pairs = target(
    command = {
      # "unwrap" command
      cdtl = ifelse(conditional == 'conditional', TRUE, FALSE)
      # generate samples
      res = lapply(standardized_sattags, function(tag) {
        list(
          # tag ID
          tag = tag$tag,
          # configuration
          response_lag = response_lag,
          window_length = window_length,
          # samples
          pairs = sample_pairs(data = data.frame(depths = tag$depths, 
                                                 depth_bins = tag$depth.bin,
                                                 diveIds = tag$diveIds),
                               times = tag$times, 
                               exposure_time = tag$exposure_time,
                               response_lag = response_lag, 
                               window_length = window_length,
                               nsamples = 1e3,
                               conditional = cdtl,
                               conditional_class = function(x) {
                                 ifelse(x['depths'] >= 800, 'Deep', 'Shallow')
                               })
        )
      })
      # save samples, and return file name
      combo = paste(res[[1]]$tag, 
                    paste('lag', response_lag, sep = ''),
                    paste('length', window_length, sep = ''),
                    sep = '_')
      f = file.path(eda_dir, paste(id_chr(), '_', combo, '.rds', 
                                   sep =''))
      saveRDS(res, file = f)
      f
    },
    dynamic = map(standardized_sattags), 
    transform = cross(response_lag = c(0, 6, 12),
                      window_length = c(1, 3, 6, 12, 18, 24, 72),
                      conditional = c('conditional', '')),
    format = 'file'
  ),
  
  # evaluate summary statistics and MC p-value for pre/post pairs
  pre_post_tests = target(
    command = {
      # get function from string
      test_fn = getFunction(test_fn_name)
      # load sampled and observed pre/post pairs
      mcout = readRDS(pre_post_pairs)[[1]]
      # initialize return
      res = mcout
      res$pairs = NULL
      res$stat = test_fn_name
      # process samples
      if(!is.null(mcout$pairs)) {
        if(length(mcout$pairs$resampled) > 0) {
          # evaluate test statistic on all pre/post samples
          res$resampled = sapply(mcout$pairs$resampled, function(sample) {
            if(is.null(sample)) {
              NA
            } else {
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
    dynamic = map(pre_post_pairs),
    transform = cross(
      pre_post_pairs, 
      test_fn_name = c('prop_surface', 'bin_counts', 'prop_downward', 
                       'prop_upward', 'prop_no_change', 'total_upward', 
                       'total_downward', 'total_absolute', 'type_seq',
                       'depth_seq', 'type_cdf', 'type_trans', 
                       'increasing_prop_downward', 'decreasing_depth_seq', 
                       'increasing_depth_seq', 'decreasing_prop_downward', 
                       'increasing_prop_no_change', 'decreasing_prop_no_change',
                       'increasing_prop_surface', 'decreasing_prop_surface', 
                       'increasing_prop_upward', 'decreasing_prop_upward', 
                       'increasing_total_absolute', 'decreasing_total_absolute',
                       'increasing_total_downward', 'decreasing_total_downward',
                       'increasing_total_upward', 'decreasing_total_upward', 
                       'increasing_type_seq', 'decreasing_type_seq')
    ),
    format = 'file'
  ),
  
  # aggregate/organize lists of pre_post files by tagId
  pre_post_files = target(
    command = list(dplyr::combine(pre_post_tests)),
    dynamic = map(pre_post_tests),
    transform = combine(pre_post_tests, .by = conditional)
  ),
  
  eda_tagsummary = target(
    command = {
      
      tag_files = pre_post_files[[1]]
      
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
      f = file.path(eda_plots_dir,
                    paste(targetTag, '_', id_chr(), '.png', sep = ''))
      ggsave(pl, filename = f)
      f
    },
    dynamic = map(pre_post_files),
    transform = map(pre_post_files),
    format = 'file'
  )

)

