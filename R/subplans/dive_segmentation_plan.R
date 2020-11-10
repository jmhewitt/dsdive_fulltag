dive_segmentation_plan = drake_plan(
  
  # location for storing sattag dive id labels
  dive_label_plot_dir = file_out(!!file.path('plots', 'tag_labels')),
  
  # location for storing dive segmentation diagnostics
  label_diagnostic_plot_dir = file_out(!!file.path('plots', 'tag_labels', 
                                                   'diagnostics')),
  
  # range of dive segmentation merge ratios to build diagnostics for
  exploratory_merge_ratios = seq(from = .005, to = .99, by = .01),
  
  # dive segmentation ratio to usse for "production"
  merge_ratio = .6,
  
  # segment tag records into dives 
  dive_labels = segment_dives(dive_label_plot_dir, label_diagnostic_plot_dir, 
                              template_bins, exploratory_merge_ratios, 
                              depth_files, merge_ratio, tag_info),
  
  # identify endpoints of segmented dives
  dive_endpoints = identify_dive_endpoints(depth_files, dive_labels, 
                                           template_bins, sattag_timestep)
  
)