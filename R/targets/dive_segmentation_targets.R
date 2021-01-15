dive_segmentation_targets = list(
  
  # range of dive segmentation merge ratios to build diagnostics for
  tar_target(exploratory_merge_ratios, seq(from = .005, to = .99, by = .01)),
  
  # dive segmentation ratio to usse for "production"
  tar_target(merge_ratio, .6),
  
  # segment tag records into dives 
  tar_target(
    name = dive_labels, 
    command = segment_dives(dive_label_plot_dir, label_diagnostic_plot_dir, 
                            template_bins, exploratory_merge_ratios, 
                            depth_files, merge_ratio, tag_info, deep_dive_depth)
  ),
  
  # identify endpoints of segmented dives
  tar_target(
    name = dive_endpoints,
    command = identify_dive_endpoints(depth_files, dive_labels, 
                                      template_bins, sattag_timestep, 
                                      deep_dive_depth)
  )
  
)