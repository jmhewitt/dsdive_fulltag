# Nimble implementation of base::expm1 via R source code, evaluating e^x - 1
expmone = nimble::nimbleFunction(
  run = function(x = double(0)) {
    
    returnType(double(0))
    
    a <- abs(x)
    
    if (a < 2.22044604925031e-16) return(x)
    
    if (a > 0.697) return(exp(x) - 1)  # negligible cancellation
      
    if (a > 1e-8) {
      y <- exp(x) - 1
    } else {
      # Taylor expansion, more accurate in this range
      y <- (x / 2 + 1) * x
    }
    
    # Newton step for solving log(1 + y) = x for y
    # WARNING: does not work for y ~ -1: bug in 1.5.0
    y <- y - (1 + y) * (log1p(y) - x)
    
    return(y)
    
  }
)

# Implement log(1 - e^x) for x < 0 via equation 7 in 
# https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
log1mexp = nimble::nimbleFunction(
  run = function(x = double(0)) {
    
    returnType(double(0))
    
    if(x < -0.693147180559945) {
      return(log1p(-exp(x)))
    }
    
    return(log(-expmone(x)))
    
  }
)

# Implement log(1 + e^x) via equation 10 in
# https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
log1pexp = nimble::nimbleFunction(
  run = function(x = double(0)) {
    
    returnType(double(0))
    
    if(x <= -37) {
      return(exp(x))
    }
    
    if(x <= 18) {
      return(log1p(exp(x)))
    }
    
    if(x <= 33.3) {
      return(x + exp(-x))
    }
    
    return(x)
    
  }
)

# implement log(c) = log(a - b) given a,b > 0, a > b, log(a), and log(b) using 
# the identity log(c) = log(a - b) = log(a) + log(1 - exp(log(b) - log(a)))
log_subtract = nimble::nimbleFunction(
  run = function(log_a = double(0), log_b = double(0)) {
    
    returnType(double(0))
    
    if(log_b == -Inf) {
      return(log_a)
    }
    
    return(log_a + log1mexp(log_b - log_a))
    
  }
)

# implement log(c) = log(a + b) given log(a) and log(b) using the identity:
#   log(c) = log(a + b) = log(b) + log( 1 + exp(log(a) - log(b)) )
log_add = nimble::nimbleFunction(
  run = function(log_a = double(0), log_b = double(0)) {
    
    returnType(double(0))
    
    if(log_a == -Inf) {
      return(log_b)
    }
    
    if(log_b == -Inf) {
      return(log_a)
    }
    
    return(log_b + log1pexp(log_a - log_b))
    
  }
)

# implements log(sum(x)) given the entrywise-logged vector log(x)
log_sum = nimble::nimbleFunction(
  run = function(log_x = double(1), n = integer(0)) {
    
    returnType(double(0))
    
    res = log_x[1]
    
    if(n == 1) {
      return(res)
    }
    
    for(i in 2:n) {
      res <- log_add(res, log_x[i])
    }
    
    return(res)
    
  }
)
