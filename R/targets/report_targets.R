report_targets = list(

  tarchetypes::tar_render(
    name = parameter_diagnostics,
    path = 'reports/templates/parameter_diagnostics.Rmd')
)
