dive_segmentation_targets = list(
  
  # range of dive segmentation merge ratios to build diagnostics for
  tar_target(exploratory_merge_ratios, seq(from = .005, to = .99, by = .01)),

  # dive segmentation ratio to usse for "production"
  tar_target(merge_ratio, .6),
  
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
  ),

  # names/labels for stages; stages are associated with alpha/beta parameters
  tar_target(
    name = stages,
    command = c('deep_descent' = 1, 'deep_forage' = 2, 'deep_ascent' = 3,
                'shallow_descent' = 4, 'shallow_ascent' = 5, 'free_surface' = 6)
  ),
  
  # map stage names to movement types (should be ordered wrt. stage labels);
  # movement types are associated with infinitesimal generator structures
  tar_target(
    name = movement_types,
    command = c('deep_descent' = 1, 'deep_forage' = 2, 'deep_ascent' = 3,
                'shallow_descent' = 1, 'shallow_ascent' = 3, 'free_surface' = 4)
  ),
  
  # factor by which to increase temporal resolution during imputation
  tar_target(imputation_factor, 5),
  
  # runs of surface obs. longer than this will be assigned to free_surface stage
  tar_target(surface_run_length, 10),
  
  # linear dive imputation
  tar_target(
    name = imputed_dive,
    command = {
      list(
      impute_observations(
        tag = raw_sattags, endpoints = dive_endpoints, 
        timestep = sattag_timestep, imputation_factor = imputation_factor, 
        template_bins = template_bins, stages = stages, 
        imputed_dive_label_plot_dir = imputed_dive_label_plot_dir,
        surface_run_length = surface_run_length
      ))
    }
    ,
    pattern = map(raw_sattags, dive_endpoints)
  )
  
)