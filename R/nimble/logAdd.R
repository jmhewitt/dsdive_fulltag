logAdd = nimble::nimbleFunction(
  run = function(log_a = double(0), log_b = double(0)) {
    # Evaluates log(a + b) when given log(a) and log(b)
    
    returnType(double(0))
    
    x <- log_a - log_b
    exp_x <- exp(x)
    
    if(exp_x == Inf) {
      res <- log_a
    } else if(exp_x == 0) {
      res <- log_b
    } else {
      res <- log_b + log(1 + exp_x)
    }
    
    return(res)
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