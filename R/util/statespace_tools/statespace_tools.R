#
# compile distribution
#

base_dir = file.path(getwd(), 'R', 'util', 'statespace_tools')
hfile = file.path(base_dir, 'statespace_tools.h')

# compile code
Rcpp::sourceCpp(file.path(base_dir,'statespace_tools.cpp'))

# retrieve compiled object file
oFile = dir(
  path = tempdir(), 
  pattern = 'statespace_tools.o', 
  full.names = TRUE, 
  recursive = TRUE
)


#
# link marginal distribution
#

# nimble wrapper to external code
marginal_ll_cpp = nimble::nimbleExternalCall(
  prototype = function(
    obs_lik_dict = double(3), obs = double(1), txmat_dict = double(3),
    txmat_seq = double(1), x0 = double(1), num_obs_states = integer(), 
    num_latent_states = integer(), nt = integer()
  ) {},
  returnType = double(),
  Cfun = 'nimLL2LayerCompressedRaw',
  headerFile = hfile,
  oFile = oFile
)

# nimble handler for external code
marginal_ll = nimble::nimbleFunction(
  run = function(
    obs_lik_dict = double(3), obs = double(1), txmat_dict = double(3),
    txmat_seq = double(1), x0 = double(1), num_obs_states = integer(0), 
    num_latent_states = integer(0), nt = integer(0)
  ) {
    returnType(double(0))
    return(
      marginal_ll_cpp(
        obs_lik_dict = obs_lik_dict, obs = obs, txmat_dict = txmat_dict, 
        txmat_seq = txmat_seq, x0 = x0, num_obs_states = num_obs_states, 
        num_latent_states = num_latent_states, nt = nt
    ))
  }
)

dstatespace = nimble::nimbleFunction(
  run = function(x = double(1), obs_lik_dict = double(3), 
                 txmat_dict = double(3), txmat_seq = double(1), x0 = double(1), 
                 num_obs_states = integer(0), num_latent_states = integer(0), 
                 nt = integer(0), log = integer(0)) {
    returnType(double(0))
    
    ll <- marginal_ll(
      obs_lik_dict = obs_lik_dict, 
      obs = x - 1, 
      txmat_dict = txmat_dict, 
      txmat_seq = txmat_seq - 1, 
      x0 = x0, 
      num_obs_states = num_obs_states,
      num_latent_states = num_latent_states,
      nt = nt
    )
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)

marginal_ll2_cpp = nimble::nimbleExternalCall(
  prototype = function(
    obs_lik_dict = double(3), obs = double(1),
    txmat_seq = double(3), x0 = double(1), num_obs_states = integer(),
    num_latent_states = integer(), nt = integer()
  ) {},
  returnType = double(),
  Cfun = 'nimLL2LayerPartialRaw',
  headerFile = hfile,
  oFile = oFile
)

marginal_ll2 = nimble::nimbleFunction(
  run = function(
    obs_lik_dict = double(3), obs = double(1),
    txmat_seq = double(3), x0 = double(1), num_obs_states = integer(0),
    num_latent_states = integer(0), nt = integer(0)
  ) {
    returnType(double(0))
    return(
      marginal_ll2_cpp(
        obs_lik_dict = obs_lik_dict, obs = obs,
        txmat_seq = txmat_seq, x0 = x0, num_obs_states = num_obs_states,
        num_latent_states = num_latent_states, nt = nt
      ))
  }
)

dstatespace2 = nimble::nimbleFunction(
  run = function(x = double(1), obs_lik_dict = double(3),
                 covariates = double(2), betas = double(3),
                 n_covariates = integer(0), x0 = double(1),
                 num_obs_states = integer(0), num_latent_states = integer(0),
                 nt = integer(0), log = integer(0)) {
    returnType(double(0))

    txmat_seq <- stageTxMats(
      betas = betas,
      covariates = covariates,
      n_timepoints = nt
    )

    ll <- marginal_ll2(
      obs_lik_dict = obs_lik_dict,
      obs = x - 1,
      txmat_seq = txmat_seq,
      x0 = x0,
      num_obs_states = num_obs_states,
      num_latent_states = num_latent_states,
      nt = nt
    )

    if(log) { return(ll) } else { return(exp(ll)) }
  }
)

nimble::registerDistributions(list(
  dstatespace = list(
    BUGSdist = paste("dstatespace(obs_lik_dict, txmat_dict, txmat_seq, x0,",
                     "num_obs_states, num_latent_states, nt)", sep = ' '),
    types = c('value = double(1)', 'obs_lik_dict = double(3)', 
              'txmat_dict = double(3)', 'txmat_seq = double(1)', 
              'x0 = double(1)', 'num_obs_states = integer(0)', 
              'num_latent_states = integer(0)', 'nt = integer(0)'),
    pqAvail = FALSE,
    discrete = TRUE
  ),
  dstatespace2 = list(
    BUGSdist = paste("dstatespace2(obs_lik_dict, covariates, betas,",
                     "n_covariates, x0, num_obs_states, num_latent_states, nt)",
                     sep = ' '),
    types = c('value = double(1)', 'obs_lik_dict = double(3)',
              'covariates = double(2)', 'n_covariates = integer(0)',
              'betas = double(3)', 'x0 = double(1)',
              'num_obs_states = integer(0)', 'num_latent_states = integer(0)',
              'nt = integer(0)'),
    pqAvail = FALSE,
    discrete = TRUE
  )
))
