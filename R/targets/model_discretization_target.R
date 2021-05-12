model_discretization_target = list(

  # discretization of parameter space
  tar_target(
    name = parameter_discretization, 
    command = list(
      # vertical speeds
      lambda = c(min_val = 0, stepsize = .05, nvals = 101),
      # directional preferences
      pi = c('ascent' = .01, 'free' = .5, 'descent' = .99)
  )),

  # transition matrices
  tar_target(
    name = transition_matrices,
    command = {

      # depth bin widths and count
      n_bins = nrow(template_bins)
      bin_widths = 2 * template_bins$halfwidth

      # range of model parameters to consider under each dive construction
      params = expand.grid(
        pi = parameter_discretization$pi,
        lambda = seq(
          from = parameter_discretization$lambda['min_val'],
          by = parameter_discretization$lambda['stepsize'],
          length.out = parameter_discretization$lambda['nvals']
        )
      )
      
      # CTMC transition matrices based on unrestricted bin movement
      stage = 6

      # compute CTMC transition matrix for this time duration
      tstep = sattag_timestep

      # condensed transition matrices for movement structures
      as.numeric(apply(params, 1, function(r) {
        expm(tstep * buildInfinitesimalGenerator(
          pi = r['pi'], lambda = r['lambda'], M = n_bins, stage = stage,
          widths = bin_widths
        ))
      }))
    }
  )

)
