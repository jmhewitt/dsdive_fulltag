multinomialLogitProbs = nimble::nimbleFunction(
  run = function(betas = double(2), x = double(1), log = logical(0)) {
    # Given covariates and coefficients, return multinomial class probabilities
    # defined via multinomial logit models, following Agresti (2002), 
    # Section 7.1.
    # 
    # Parameters:
    #   betas - p x (J-1) matrix of coefficients for multinomial logits.  
    #    The ith column holds coefficients for logit(pi_i/pi_J) where J is 
    #    the reference class and i = 1, ..., J-1.
    #   x -  p-dimensional vector of covariates to be used with coefficients.
    #   log - if TRUE, returns an approximation to the log probabilities
    
    returnType(double(1))
    
    # eqn. 7.1, with alpha_j assumed to live in betas as an intercept term
    logits <- c(t(betas) %*% x, 0)
    exp_logits <- exp(logits)
    
    
    # eqn. 7.2
    norm_const <- sum(exp_logits)
    if(log) {
      res <- logits - log(norm_const)
    } else {
      res <- exp_logits / norm_const
    }
    
    return(res)
  }
)