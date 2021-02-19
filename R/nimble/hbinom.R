hbinom = nimble::nimbleFunction(
  run = function(x = double(0), size = double(0), prob = double(0), 
                 log = logical(0, default = 0)) {
    # hazard function for binomial distribution.
    #
    # Parameters:
    #  x - sequence of states
    
    returnType(double(0))
    
    # evaluate log-hazard function
    if(x > size) {
      h<- 0 
    } else {
      h <- dbinom(x = x, size = size, prob = prob, log = TRUE) -
        pbinom(q = x, size = size, prob = prob, lower.tail = FALSE, log.p = TRUE)
    }
    
    if(log) { return(h) } else { return(exp(h))}
  }
)
