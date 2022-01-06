dive_segmentation_targets = list(
  
  # range of dive segmentation merge ratios to build diagnostics for
  tar_target(exploratory_merge_ratios, seq(from = .005, to = .99, by = .01)),

  # dive segmentation ratio to usse for "production"
  tar_target(merge_ratio, .6),
  
  # location for storing sattag dive id labels
  tar_target(
    name = dive_label_plot_dir, 
    command = {
      f = file.path('output', 'tag_labels')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # location for storing dive segmentation diagnostics
  tar_target(
    name = label_diagnostic_plot_dir, 
    command = {
      f = file.path(dive_label_plot_dir, 'diagnostics')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # segment tag records into dives 
  tar_target(
    name = dive_labels, 
    command = segment_dives(
      dive_label_plot_dir = dive_label_plot_dir, 
      label_diagnostic_plot_dir = label_diagnostic_plot_dir, 
      template_bins = template_bins, 
      exploratory_merge_ratios = exploratory_merge_ratios, 
      depth_files = depth_files, 
      merge_ratio = merge_ratio, 
      tag_info = tag_info, deep_threshold = deep_dive_depth, 
      timestep = sattag_timestep
    )
  ),
  
  # identify endpoints of segmented dives
  tar_target(
    name = dive_endpoints,
    command = identify_dive_endpoints(depth_files = depth_files, 
                                      dive_labels = dive_labels, 
                                      template_bins = template_bins, 
                                      sattag_timestep = sattag_timestep, 
                                      deep_threshold = deep_dive_depth)
  )
  
)