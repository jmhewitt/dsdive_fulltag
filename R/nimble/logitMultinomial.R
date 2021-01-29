multinomialLogitProbs = nimble::nimbleFunction(
  run = function(betas = double(2), x = double(1)) {
    # Given covariates and coefficients, return multinomial class probabilities
    # defined via multinomial logit models, following Agresti (2002), 
    # Section 7.1.
    # 
    # Parameters:
    #   betas - p x (J-1) matrix of coefficients for multinomial logits.  
    #    The ith column holds coefficients for logit(pi_i/pi_J) where J is 
    #    the reference class and i = 1, ..., J-1.
    #   x -  p-dimensional vector of covariates to be used with coefficients.
    
    returnType(double(1))
    
    # eqn. 7.1, with alpha_j assumed to live in betas as an intercept term
    exp_logits <- exp(c(t(betas) %*% x, 1))
    
    # eqn. 7.2
    prob <- exp_logits / sum(exp_logits)
    
    return(prob)
  }
)