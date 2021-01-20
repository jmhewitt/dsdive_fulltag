lookupProb = nimble::nimbleFunction(
  run = function(movement_type = integer(0), pi_ind = integer(0), 
                 lambda_ind = integer(0), i = integer(0), j = integer(0), 
                 n_pi = integer(1), n_lambda = integer(1), n_bins = integer(0), 
                 tmats = double(1)) {
    # Parameters:
    #  movement_type - the type of movement (i.e., ascent, forage, descent) 
    #    the probability lookup is for
    #  pi_ind - the index of the pre-specified descent preference value
    #  lambda_ind - the index of the pre-specified speed value
    #  i - the depth bin being transitioned from
    #  j - the depth bin being transitioned to
    #  n_pi - the number of pre-specified descent preference values
    #  n_lambda - the number of pre-specified speed values
    #  n_bins - the number of depth bins
    #  tmats - the transition matrices, stored in flat format
    
    returnType(double(0))
    
    # elements w/in a single transition matrix
    Ksq <- n_bins^2
    
    # index w/in target matrix (stored in col. major format)
    index <- n_bins * (j-1) + i
    
    # offset wrt. movement type
    if(movement_type > 1) {
      index <- index + Ksq * n_pi[1] * n_lambda[1]
    }
    if(movement_type > 2) {
      index <- index + Ksq * n_pi[2] * n_lambda[2]
    }
    
    # return after offsetting wrt. grouping vars
    p <- tmats[
      index + Ksq * (n_pi[movement_type] * (lambda_ind - 1) + pi_ind - 1)
      ]
    
    return(p)
  }
)


# Demo to validate the methods work as intended
# 
# mt = 1
# pi_ind = 26
# lambda_ind = 1
# i = 3
# j = 2
# 
# x = lookupProb(movement_type = mt, pi_ind = pi_ind, 
#                lambda_ind = lambda_ind, i = i, j = j, 
#                n_pi = parameter_discretization$pi[, 'nvals'], 
#                n_lambda = parameter_discretization$lambda[, 'nvals'], 
#                n_bins = n_bins, tmats = tmats)
# 
# xm = lookupTmat(movement_type = mt, pi_ind = pi_ind, 
#                 lambda_ind = lambda_ind, 
#                 n_pi = parameter_discretization$pi[, 'nvals'], 
#                 n_lambda = parameter_discretization$lambda[, 'nvals'], 
#                 n_bins = n_bins, tmats = tmats)
# 
# validateMat = expm(tstep * buildInfinitesimalGenerator(
#   pi = seq(from = parameter_discretization$pi[mt, 'min_val'], 
#            by = parameter_discretization$pi[mt, 'stepsize'], 
#            length.out = parameter_discretization$pi[mt, 'nvals'])[pi_ind], 
#   lambda = seq(from = parameter_discretization$lambda[mt, 'min_val'], 
#                by = parameter_discretization$lambda[mt, 'stepsize'], 
#                length.out = parameter_discretization$lambda[mt, 'nvals'])[lambda_ind], 
#   M = n_bins, 
#   stage = mt, 
#   widths = bin_widths
# ))
# 
# identical(x, validateMat[i,j])
# identical(xm, validateMat)