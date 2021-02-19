hnbinom = nimble::nimbleFunction(
  run = function(x = double(0), size = double(0), prob = double(0), 
                 log = logical(0, default = 0), shift = double(0)) {
    # hazard function for shifted negative-binomial distribution.  evaluates 
    # hazard function at h(x - shift).
    #
    # Parameters:
    #  x - sequence of states
    
    returnType(double(0))
    
    # evaluate log-hazard function
    h <- dnbinom(x = x - shift, size = size, prob = prob, log = TRUE) - 
      pnbinom(q = x - shift, size = size, prob = prob, lower.tail = FALSE, 
              log.p = TRUE)
    
    if(log) { return(h) } else { return(exp(h))}
  }
)
