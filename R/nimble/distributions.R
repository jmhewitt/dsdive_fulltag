cpp.dir = file.path('scripts', 'nimble')

#
# support for CTMC model across depth bins
#

# nimbleList for storing entries of a tridiagonal matrix
tridiagonalEntries = nimbleList(
  diag = double(1),
  dsuper = double(1),
  dsub = double(1)
)

# build components for infinitesimal generator matrix for depth bin transitions
buildInfinitesimalGeneratorEntries = nimbleFunction(
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
buildInfinitesimalGenerator = nimbleFunction(
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


#
# C++ support for matrix exponentials
#

# compile C++ implementation of LAPACK-assisted generator decomposition
system(paste(
  'R CMD COMPILE',
  paste('CXXFLAGS=', '"',
        paste('-I', 
              file.path(find.package(c('RcppEigen', 'Rcpp')), 'include'), 
              sep = '', collapse = ' '), 
        '"', sep = ''),
  file.path(cpp.dir, 'exp_symtools.cpp'),
  sep = ' '
))

# link NIMBLE to external C++ (decomposition)
decomposeGeneratorCpp = nimbleExternalCall(
  prototype = function(diag = double(1), dsuper = double(1), dsub = double(1), 
                       N = integer(0), expm = double(2), evecs = double(2), 
                       d = double(1), dInv = double(1), delta = double(0),
                       t = double(0)){},  
  returnType = void(), 
  Cfun = 'expm_cpp', 
  headerFile = file.path(getwd(), cpp.dir, 'exp_symtools.h'), 
  oFile = file.path(getwd(), cpp.dir, 'exp_symtools.o')
)

# nimbleList for eigendecomposition of matrix exponential
eigenExpm = nimbleList(
  expm = double(2),
  evecs = double(2),
  evals = double(1),
  d = double(1),
  dInv = double(1),
  t = double(0),
  N = integer(0)
)

# construct eigendecomposition of matrix exponential
decomposeGenerator = nimbleFunction(
  run = function(Aentries = tridiagonalEntries(), delta = double(0), 
                 t = double(0), N = integer(0)) {

    returnType(eigenExpm())
    
    expm <- matrix(0, nrow = N, ncol = N)
    evecs <- matrix(0, nrow = N, ncol = N)
    d <- numeric(0, length = N)
    dInv <- numeric(0, length = N)
    
    evals <- numeric(length = N, init = FALSE)
    evals[1:N] <- Aentries$diag[1:N]

    decomposeGeneratorCpp(diag = evals, dsuper = Aentries$dsuper, 
                          dsub = Aentries$dsub, N = N, expm = expm, 
                          evecs = evecs, d = d, dInv = dInv, delta = delta, 
                          t = t)
    
    res <- eigenExpm$new(expm = expm, evecs = evecs, evals = evals, d = d,
                         dInv = dInv, t = t, N = N)
    
    return(res)
  }
)

# link NIMBLE to external C++ (multiplication)
expmAtvCpp = nimbleExternalCall(
  prototype = function(evecs = double(2), evals = double(1), M = integer(0), 
                       v = double(1), d = double(1), dInv = double(1), 
                       t = double(0), x = double(1), 
                       preMultiply = logical(0)){},
  returnType = void(), 
  Cfun = 'expmAtv_cpp', 
  headerFile = file.path(getwd(), cpp.dir, 'exp_symtools.h'), 
  oFile = file.path(getwd(), cpp.dir, 'exp_symtools.o')
)

# use generator matrix decomposition to compute x = exp(At)*v
expmAtv = nimbleFunction(
  run = function(expm = double(2), evecs = double(2), evals = double(1), 
                 d = double(1), dInv = double(1), tstep = double(0), 
                 N = integer(0), v = double(1), t = double(0), 
                 preMultiply = logical(0)) {
    
    returnType(double(1))
    
    x <- numeric(0, length = N)
    
    if(t == tstep) {
      if(preMultiply == TRUE) {
        x[1:N] <- t(v[1:N]) %*% expm[1:N,1:N]
      } else {
        x[1:N] <- expm[1:N,1:N] %*% v[1:N]
      }
    } else {
      expmAtvCpp(evecs = evecs, evals = evals, M = N, v = v, d = d, dInv = dInv, 
                 t = t, x = x, preMultiply = preMultiply)
    }
    
    return(x)
  }
)

buildAndDecomposeGenerator = nimbleFunction(
  run = function(pi = double(0), lambda = double(0), M = integer(0),
                 stage = integer(0), widths = double(1), delta = double(0), 
                 t = double(0)) {
    
    returnType(double(2))
    
    res <- matrix(nrow = 2 * M + 3, ncol = M, init = FALSE)
    
    # build generator entries
    entries <- buildInfinitesimalGeneratorEntries(
      pi = pi, lambda = lambda, M = M, stage = stage, widths = widths
    )
    
    # decompose generator
    decomposition <- decomposeGenerator(Aentries = entries, delta = delta, 
                                        t = t, N = M)
    
    # assign output
    res[1:M,1:M] <- decomposition$expm[1:M,1:M]
    res[(M+1):(2*M),1:M] <- decomposition$evecs[1:M,1:M]
    res[2*M+1, 1:M] <- decomposition$evals[1:M]
    res[2*M+2,1:M] <- decomposition$d[1:M]
    res[2*M+3,1:M] <- decomposition$dInv[1:M]
    
    return(res)
  }
)


#
# prior densities
#

# density for logit(X) where X ~ Beta(shape1, shape2)
dlogitBeta = nimbleFunction(
  run = function(x = double(0), shape1 = double(0), shape2 = double(0), 
                 log = logical(0, default = 0)) {
    
    returnType(double(0))
    
    res <- dbeta(x = ilogit(x), shape1 = shape1, shape2 = shape2, log = TRUE) + 
      x - 2 * log(exp(x) + 1)
    
    if(log) { return(res) } else { return(exp(res)) }
  }
)

# random generation for logit(X) where X ~ Beta(shape1, shape2)
rlogitBeta = nimbleFunction(
  run = function(n = integer(0), shape1 = double(0), shape2 = double(0)) {
    returnType(double(0))
    res <- rbeta(n = 1, shape1 = shape1, shape2 = shape2)
    return(logit(res))
  }
)

registerDistributions(list(
  dlogitBeta = list(
    BUGSdist = 'dlogitBeta(shape1, shape2)',
    Rdist = 'dlogitBeta(shape1, shape2)',
    types = c('shape1 = double(0)', 'shape2 = double(0)')
  )
))

# density for log(X) where X ~ Gamma(shape, rate)
dlogGamma = nimbleFunction(
  run = function(x = double(0), shape = double(0), rate = double(0), 
                 log = logical(0, default = 0)) {
    
    returnType(double(0))
    
    res <- dgamma(x = exp(x), shape = shape, rate = rate, log = TRUE) + x
    
    if(log) { return(res) } else { return(exp(res)) }
  }
)

#
# likelihood function
#

# likelihood for an entire deep dive
ddive = nimbleFunction(
  run = function(x = double(1), times = double(1), expm = double(3), 
                 N = integer(0),
                 evecs = double(3), evals = double(2), d = double(2), 
                 dInv = double(2), tstep = double(0), include = logical(0),
                 M = integer(0), T = double(1), log = logical(0, default = 0)) {
    # Parameters:
    #   x - sequence of observed depth bins 
    #   times - observation times
    #   N - number of observations
    #   M - number of depth bins
    #   T - vector of dive-specific random effects;  c(T0, T1, T2, T3) s.t. 
    #       T0 is dive start time, T1 and T2 are respecitively end of stage 1/2,
    #       and T3 is dive end time.  note that the indexing implies T[1] is 
    #       start of dive, T[2] is start of stage 2, T[3] is start of stage 3, 
    #       and T[4] is start of surface period (i.e., end of dive).
    #   log - If TRUE, then log-likelihood is returned
    
    returnType(double(0))
  
    if(include == FALSE) {
      ll <- 0
    } else {
      # density is not defined for mis-ordered dive-specific random effects
      if(T[1] >= T[2]) { return(-Inf) }
      if(T[2] >= T[3]) { return(-Inf) }
      if(T[3] >= T[4]) { return(-Inf) }

      # initialize log-likelihood
      ll <- 0

      # initialize Markov transitions from beginning of dive
      x_prev <- 1
      time_prev <- T[1]
      stage_prev <- 1

      # aggregate likelihood until end of dive
      more_observations <- TRUE
      i <- 1
      while(more_observations) {

        if(i <= N) {
          # ll contribution from next observation
          time_next <- times[i]
          x_next <- x[i]
          i <- i+1
        } else {
          # ll contribution from dive end
          time_next <- T[4]
          x_next <- 1
          more_observations <- FALSE
        }

        # only process in-sync observations (e.g., times within dive window)
        if(time_next > time_prev) {

          # extract current transition, with adjustment for end of dive
          if(time_next >= T[4]) {
            x_loc <- 1
            time_loc <- T[4]
            more_observations <- FALSE
          } else {
            x_loc <- x_next
            time_loc <- time_next
          }

          # determine dive stage at end of transition
          if(time_loc >= T[3]) {
            stage_loc <- 3
          } else if(time_loc >= T[2]) {
            stage_loc <- 2
          } else {
            stage_loc <- 1
          }

          # initialize transition distribution
          u <- numeric(0, length = M)
          u[x_prev] <- 1

          # diffuse transition mass across stages
          time_stage <- time_prev
          for(s in stage_prev:stage_loc) {
            # time spent in stage
            d_t <- min(T[s+1], time_loc) - time_stage
            # update time-in-stage marker
            time_stage <- T[s+1]
            # diffuse mass across stage s
            u[1:M] <- expmAtv(expm = expm[s,1:M,1:M], evecs = evecs[s,1:M,1:M],
                              evals = evals[s,1:M], d = d[s,1:M],
                              dInv = dInv[s,1:M], tstep = tstep, N = M,
                              v = u[1:M], t = d_t, preMultiply = TRUE)[1:M]
            mass <- sum(u[1:M])
            for(k in 1:M) {
              u[k] <- u[k] / mass
            }
          }

          # aggregate likelihood
          ll <- ll + log(u[x_loc])

          # update last state, for computing Markov probabilities
          x_prev <- x_loc
          time_prev <- time_loc
          stage_prev <- stage_loc
        }
      }
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)


# likelihood for an entire shallow dive
ddiveShallow = nimbleFunction(
  run = function(x = double(1), times = double(1), expm = double(3), 
                 N = integer(0),
                 evecs = double(3), evals = double(2), d = double(2), 
                 dInv = double(2), tstep = double(0), include = logical(0),
                 M = integer(0), T = double(1), log = logical(0, default = 0)) {
    # Parameters:
    #   x - sequence of observed depth bins 
    #   times - observation times
    #   N - number of observations
    #   M - number of depth bins
    #   T - vector of dive-specific random effects;  c(T0, T1, T2, T3) s.t. 
    #       T0 is dive start time, T1 is the end of stage 1, T2 is not used,
    #       and T3 is dive end time.  note that the indexing implies T[1] is 
    #       start of dive, T[2] is start of stage 2, 
    #       and T[4] is start of surface period (i.e., end of dive).
    #   log - If TRUE, then log-likelihood is returned
    
    returnType(double(0))
    
    if(include == FALSE) {
      ll <- 0
    } else {
      # density is not defined for mis-ordered dive-specific random effects
      if(T[1] >= T[2]) { return(-Inf) }
      if(T[2] >= T[4]) { return(-Inf) }
      
      # simplify for later use
      T[3] <- T[4]
      
      # initialize log-likelihood
      ll <- 0
      
      # initialize Markov transitions from beginning of dive
      x_prev <- 1
      time_prev <- T[1]
      stage_prev <- 1
      
      # aggregate likelihood until end of dive
      more_observations <- TRUE
      i <- 1
      while(more_observations) {
        
        if(i <= N) {
          # ll contribution from next observation
          time_next <- times[i]
          x_next <- x[i]
          i <- i+1
        } else {
          # ll contribution from dive end
          time_next <- T[4]
          x_next <- 1
          more_observations <- FALSE
        }
        
        # only process in-sync observations (e.g., times within dive window)
        if(time_next > time_prev) {
          
          # extract current transition, with adjustment for end of dive
          if(time_next >= T[4]) {
            x_loc <- 1
            time_loc <- T[4]
            more_observations <- FALSE
          } else {
            x_loc <- x_next
            time_loc <- time_next
          }
          
          # determine dive stage at end of transition
          if(time_loc >= T[2]) {
            stage_loc <- 2
          } else {
            stage_loc <- 1
          }
          
          # initialize transition distribution
          u <- numeric(0, length = M)
          u[x_prev] <- 1
          
          # diffuse transition mass across stages
          time_stage <- time_prev
          for(s in stage_prev:stage_loc) {
            # time spent in stage
            d_t <- min(T[s+1], time_loc) - time_stage
            # update time-in-stage marker
            time_stage <- T[s+1]
            # diffuse mass across stage s
            u[1:M] <- expmAtv(expm = expm[s,1:M,1:M], evecs = evecs[s,1:M,1:M], 
                              evals = evals[s,1:M], d = d[s,1:M], 
                              dInv = dInv[s,1:M], tstep = tstep, N = M, 
                              v = u[1:M], t = d_t, preMultiply = TRUE)[1:M]
            mass <- sum(u[1:M])
            for(k in 1:M) {
              u[k] <- u[k] / mass
            }
          }
          
          # aggregate likelihood
          ll <- ll + log(u[x_loc])
          
          # update last state, for computing Markov probabilities
          x_prev <- x_loc
          time_prev <- time_loc
          stage_prev <- stage_loc
        }
      }
    }
      
      if(log) { return(ll) } else { return(exp(ll)) }
    }
)


ddive = nimbleFunction(
  run = function(x = double(1), times = double(1), expm = double(3), 
                 N = integer(0),
                 evecs = double(3), evals = double(2), d = double(2), 
                 dInv = double(2), tstep = double(0), include = logical(0),
                 M = integer(0), T = double(1), log = logical(0, default = 0)) {
    # Parameters:
    #   x - sequence of observed depth bins 
    #   times - observation times
    #   N - number of observations
    #   M - number of depth bins
    #   T - vector of dive-specific random effects;  c(T0, T1, T2, T3) s.t. 
    #       T0 is dive start time, T1 and T2 are respecitively end of stage 1/2,
    #       and T3 is dive end time.  note that the indexing implies T[1] is 
    #       start of dive, T[2] is start of stage 2, T[3] is start of stage 3, 
    #       and T[4] is start of surface period (i.e., end of dive).
    #   log - If TRUE, then log-likelihood is returned
    
    returnType(double(0))
    
    if(include == FALSE) {
      ll <- 0
    } else {
      # density is not defined for mis-ordered dive-specific random effects
      if(T[1] >= T[2]) { return(-Inf) }
      if(T[2] >= T[3]) { return(-Inf) }
      if(T[3] >= T[4]) { return(-Inf) }
      
      # initialize log-likelihood
      ll <- 0
      
      # initialize Markov transitions from beginning of dive
      x_prev <- 1
      time_prev <- T[1]
      stage_prev <- 1
      
      # aggregate likelihood until end of dive
      more_observations <- TRUE
      i <- 1
      while(more_observations) {
        
        if(i <= N) {
          # ll contribution from next observation
          time_next <- times[i]
          x_next <- x[i]
          i <- i+1
        } else {
          # ll contribution from dive end
          time_next <- T[4]
          x_next <- 1
          more_observations <- FALSE
        }
        
        # only process in-sync observations (e.g., times within dive window)
        if(time_next > time_prev) {
          
          # extract current transition, with adjustment for end of dive
          if(time_next >= T[4]) {
            x_loc <- 1
            time_loc <- T[4]
            more_observations <- FALSE
          } else {
            x_loc <- x_next
            time_loc <- time_next
          }
          
          # determine dive stage at end of transition
          if(time_loc >= T[3]) {
            stage_loc <- 3
          } else if(time_loc >= T[2]) {
            stage_loc <- 2
          } else {
            stage_loc <- 1
          }
          
          # initialize transition distribution
          u <- numeric(0, length = M)
          u[x_prev] <- 1
          
          # diffuse transition mass across stages
          time_stage <- time_prev
          for(s in stage_prev:stage_loc) {
            # time spent in stage
            d_t <- min(T[s+1], time_loc) - time_stage
            # update time-in-stage marker
            time_stage <- T[s+1]
            # diffuse mass across stage s
            u[1:M] <- expmAtv(expm = expm[s,1:M,1:M], evecs = evecs[s,1:M,1:M],
                              evals = evals[s,1:M], d = d[s,1:M],
                              dInv = dInv[s,1:M], tstep = tstep, N = M,
                              v = u[1:M], t = d_t, preMultiply = TRUE)[1:M]
            mass <- sum(u[1:M])
            for(k in 1:M) {
              u[k] <- u[k] / mass
            }
          }
          
          # aggregate likelihood
          ll <- ll + log(u[x_loc])
          
          # update last state, for computing Markov probabilities
          x_prev <- x_loc
          time_prev <- time_loc
          stage_prev <- stage_loc
        }
      }
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)


function(x = double(1), times = double(1), expm = double(3), 
         N = integer(0),
         evecs = double(3), evals = double(2), d = double(2), 
         dInv = double(2), tstep = double(0), include = logical(0),
         M = integer(0), T = double(1), log = logical(0, default = 0))
  
  
# dummy simulation function for ddive; needed to make drake happy with nimble
rdive = nimbleFunction(
  run = function (n = integer(0), times = double(1), expm = double(3), 
                  N = integer(0), evecs = double(3), evals = double(2), 
                  d = double(2), dInv = double(2), tstep = double(0), 
                  include = logical(0), M = integer(0), T = double(1)) {
    returnType(double(1))
  x <- nimNumeric(length = n)
  return(x)
})

registerDistributions(list(
  ddive = list(
    BUGSdist = 
      'ddive(times, expm, N, evecs, evals, d, dInv, tstep, include, M, T)',
    types = c('value = double(1)', 'times = double(1)', 'expm = double(3)', 
              'N = integer(0)', 'evecs = double(3)', 'evals = double(2)', 
              'd = double(2)', 'dInv = double(2)', 'tstep = double(0)', 
              'include = logical(0)', 'M = integer(0)', 'T = double(1)')
  )
))
