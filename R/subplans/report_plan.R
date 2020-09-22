report_plan = drake_plan(
  
  posterior_report = target(
    rmarkdown::render(
      knitr_in(!!file.path('reports', 'templates', 
                           'posterior_diagnostics.Rmd')),
      output_file = 'posterior_diagnostics.html',
      output_dir = 'reports',
      quiet = FALSE,
      params = list(samples = mcmc_samples,
                    nim_pkg_tgt = 'nim_pkg',
                    nburn = nburn)
    )
  ),
  
  # posterior_report_stagelearning = target(
  #   rmarkdown::render(
  #     knitr_in(!!file.path('reports', 'templates', 
  #                          'posterior_diagnostics.Rmd')),
  #     output_file = 'posterior_diagnostics_stagelearning.html',
  #     output_dir = 'reports',
  #     quiet = FALSE,
  #     params = list(samples = mcmc_samples_stagelearning,
  #                   nim_pkg_tgt = 'nim_pkg',
  #                   nburn = nburn)
  #   )
  # ),
  
  posterior_durations_stagelearning = 
    target(
      rmarkdown::render(
        knitr_in(!!file.path('reports', 'templates', 
                             'posterior_stagelearning.Rmd')),
        output_file = 'posterior_durations_stagelearning.html',
        output_dir = 'reports',
        quiet = FALSE,
        params = list(samples = mcmc_samples_stagelearning,
                      nim_pkg_tgt = 'nim_pkg',
                      nburn = nburn)
      ), 
      trigger = trigger(condition = FALSE)
    ),
  
  exposure_eda = 
    target(
      rmarkdown::render(
        knitr_in(!!file.path('reports', 'exposure_eda.Rmd')),
        output_file = paste(.id_chr, '.html', sep = ''),
        output_dir = file.path('reports', 'exposure_eda'),
        quiet = FALSE,
        params = list(window_length = window_length, 
                      response_lag = response_lag)
      ), 
      trigger = trigger(condition = TRUE), 
      transform = cross(window_length = c(1, 3, 6, 12, 18, 24, 72),
                        response_lag = c(0, 6, 12))
    )
  
)

