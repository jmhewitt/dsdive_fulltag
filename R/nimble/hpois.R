hpois = nimble::nimbleFunction(
  run = function(x = double(0), lambda = double(0), 
                 log = logical(0, default = 0)) {
    # hazard function for poisson distribution.
    #
    # Parameters:
    #  x - sequence of states
    
    returnType(double(0))
    
    # evaluate log-hazard function
    h <- dpois(x = x, lambda = lambda, log = TRUE) -
      ppois(q = x, lambda = lambda, lower.tail = FALSE, log.p = TRUE)
    
    if(log) { return(h) } else { return(exp(h))}
  }
)
