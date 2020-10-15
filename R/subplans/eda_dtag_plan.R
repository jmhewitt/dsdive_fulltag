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
    load_raw_dtag(dtag_files = dtag_files, cee_starts = cee_starts),
    dynamic = map(dtag_files)
  )
  

)

