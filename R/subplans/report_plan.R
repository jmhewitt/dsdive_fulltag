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
      trigger = trigger(condition = TRUE)
    )
  
)

