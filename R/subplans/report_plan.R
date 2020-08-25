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
  )
  
)

