#
# support for CTMC model across depth bins
#

# nimbleList for storing entries of a tridiagonal matrix
tridiagonalEntries = nimble::nimbleList(
  diag = double(1),
  dsuper = double(1),
  dsub = double(1)
)

# build components for infinitesimal generator matrix for depth bin transitions
buildInfinitesimalGeneratorEntries = nimble::nimbleFunction(
  run = function(pi = double(0), lambda = double(0), M = integer(0),
                 stage = integer(0), widths = double(1)) {
    # Parameters:
    #  pi - probability of making a downward descent
    #  lambda - speed of moving through bin
    #  M - total number of depth bins
    #  stage - build generator matrix for "stage"
    #  widths - vector of depth bin widths

    returnType(tridiagonalEntries())

    entries <- tridiagonalEntries$new(diag = numeric(M, init = FALSE),
                                      dsub = numeric(M-1, init = FALSE),
                                      dsuper = numeric(M-1, init = FALSE))

    # only allow (downward) transitions from shallowest bin during stages 1 & 2
    if(stage == 1 | stage == 2) {
      rate <- lambda / widths[1]
      entries$diag[1] <- -rate
      entries$dsuper[1] <- rate
    }

    # only allow surfacing from second bin to occur in stage 3
    rate <- lambda / widths[2]
    if(stage == 1 | stage == 2) {
      entries$diag[2] <- -rate
      entries$dsuper[2] <- rate
      entries$dsub[1] <- 0
    } else {
      entries$diag[2] <- -rate
      entries$dsuper[2] <- rate * pi
      entries$dsub[1] <- rate * (1-pi)
    }

    # intermediate bins may always transition up or down
    for(i in 3:(M-1)) {
      rate <- lambda / widths[i]
      entries$diag[i] <- -rate
      entries$dsuper[i] <- rate * pi
      entries$dsub[i-1] <- rate * (1-pi)
    }

    # deepest bin can only transition to next shallowest depth
    rate <- lambda / widths[M]
    entries$diag[M] <- -rate
    entries$dsub[M-1] <- rate

    return(entries)
  }
)

# build infinitesimal generator matrix for depth bin transitions
buildInfinitesimalGenerator = nimble::nimbleFunction(
  run = function(pi = double(0), lambda = double(0), M = integer(0),
                 stage = integer(0), widths = double(1)) {
    # Parameters:
    #  pi - probability of making a downward descent
    #  lambda - speed of moving through bin
    #  M - total number of depth bins
    #  stage - build generator matrix for "stage"
    #  widths - vector of depth bin widths

    returnType(double(2))

    # compute matrix's tridiagonal entries
    entries <- buildInfinitesimalGeneratorEntries(pi = pi, lambda = lambda,
                                                  M = M, stage = stage,
                                                  widths = widths)

    # initialize matrix
    A <- matrix(0, nrow = M, ncol = M)

    #
    # populate matrix
    #

    for(i in 1:(M-1)) {
      A[i,i] <- entries$diag[i]
      A[i,i+1] <- entries$dsuper[i]
      A[i+1,i] <- entries$dsub[i]
    }

    A[M,M] <- entries$diag[M]

    return(A)
  }
)