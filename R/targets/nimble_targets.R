nimble_targets = list(
  
  tar_target(
    name = nim_pkg, 
    command = flatten_tags(
      tag_list = imputed_dive, 
      transition_matrices = transition_matrices, 
      n_bins = nrow(template_bins), 
      movement_types = movement_types, 
      pi_discretization = parameter_discretization$pi, 
      lambda_discretization = parameter_discretization$lambda
    )
  )
  
)
