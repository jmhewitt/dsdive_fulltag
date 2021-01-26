model_discretization_target = list(

  # discretization of parameter space
  tar_target(parameter_discretization, list(

    # vertical speeds, by movement structure (i.e., stage)
    lambda = matrix(
      rep(c(0, .05, 101), 3),
      ncol = 3, byrow = TRUE,
      dimnames = list(
        rownames = c('descent', 'foraging', 'ascent'),
        colnames = c('min_val', 'stepsize', 'nvals')
      )
    ),

    # directional preferences, by movement structure (i.e., stage)
    pi = matrix(
      c(0, .01, 51,
        0.5, 0, 1,
        0.5, .01, 51),
      ncol = 3, byrow = TRUE,
      dimnames = list(
        rownames = c('descent', 'foraging', 'ascent'),
        colnames = c('min_val', 'stepsize', 'nvals')
      )
    )
  )),

  # transition matrices
  tar_target(
    name = transition_matrices,
    command = {

      # depth bin widths and count
      n_bins = nrow(template_bins)
      bin_widths = 2 * template_bins$halfwidth

      # range of model parameters to consider under each dive construction
      movement_types = c('descent', 'foraging', 'ascent')
      params = lapply(movement_types, function(mv_type) {
        expand.grid(
          pi = seq(
            from = parameter_discretization$pi[mv_type, 'min_val'],
            by = parameter_discretization$pi[mv_type, 'stepsize'],
            length.out = parameter_discretization$pi[mv_type, 'nvals']
          ),
          lambda = seq(
            from = parameter_discretization$lambda[mv_type, 'min_val'],
            by = parameter_discretization$lambda[mv_type, 'stepsize'],
            length.out = parameter_discretization$lambda[mv_type, 'nvals']
          )
        )
      })
      names(params) = movement_types

      # compute CTMC transition matrix for this time duration
      tstep = sattag_timestep

      # condensed transition matrices for movement structures
      as.numeric(do.call(cbind, mapply(function(params, stage) {
        apply(params, 1, function(r) {
          expm(tstep * buildInfinitesimalGenerator(
            pi = r['pi'], lambda = r['lambda'], M = n_bins, stage = stage,
            widths = bin_widths
          ))
        })
      }, params, 1:3, SIMPLIFY = FALSE)))
    }
  )

)
