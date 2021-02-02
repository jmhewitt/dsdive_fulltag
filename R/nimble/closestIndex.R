closestIndex = nimble::nimbleFunction(
  run = function(value = double(0), minval = double(0), stepsize = double(0),
                 nvals = double(0)) {
    # Parameters:
    #  value - continuous-space parameter to map to discrete set
    #  minval - minimum value in discrete set
    #  stepsize - spacing between values in discrete set
    #  nvals - number of values in set
    
    returnType(integer(0))
    
    # determine "optimal" index
    ind <- round((value - minval) / stepsize) + 1

    # constrain optimal index to be within the discrete set
    if(ind < 0) {
      ind <- 1
    } else if(ind > nvals) {
      ind <- nvals
    }
    
    return(ind)
  }
)
