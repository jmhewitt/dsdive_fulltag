lognormal.params = function(mean, var) {
  # Parameters:
  #  mean - E(X) for X ~ log-Normal(mu, sigmasq)
  #  var - Var(X) for the same
  # 
  # Output:
  #  Parameters c(mu, sigmasq) for a log-Normal distribution that matches given
  #  mean and variance
  
  sigmasq = log(var / mean^2 + 1)
  c(mu = log(mean) - sigmasq/2, sigmasq = sigmasq)
}
