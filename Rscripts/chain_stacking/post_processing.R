library(loo)
library(coda)
library(stringr)

source(file.path('R', 'util', 'log_add.R'))

c_log_subtract = nimble::compileNimble(log_subtract)

ncores = parallel::detectCores()/2

# identify files
f = dir(
  path = file.path('output', 'mcmc', 'chain_weights'), 
  pattern = 'model_ll_samples.rds', 
  full.names = TRUE, 
  recursive = TRUE
)

# load likelihood samples
ll_list = lapply(f, function(f) readRDS(f)[,-1])

# transform to array
ll_array = array(
  data = NA, 
  dim = c(nrow(ll_list[[1]]), length(ll_list), ncol(ll_list[[1]]))
)
for(i in 1:length(ll_list)) {
  ll_array[,i,] = ll_list[[i]]
}

r_eff = relative_eff(x = ll_array, cores = ncores)
  
loo_list = lapply(ll_list, function(x) {
  r = loo(x, r_eff = r_eff, cores = ncores)
})

lpd_point <- do.call(
  cbind, lapply(loo_list, function(obj) obj$pointwise[, "elpd_loo"])
)

stacking_weights_local = function(
    lpd_point, optim_method = "BFGS", optim_control = list()
) {
  stopifnot(is.matrix(lpd_point))
  N <- nrow(lpd_point)
  K <- ncol(lpd_point)
  if (K < 2) {
    stop("At least two models are required for stacking weights.")
  }
  exp_lpd_point <- exp(lpd_point)
  negative_log_score_loo <- function(w) {
    stopifnot(length(w) == K - 1)
    w_full <- log(c(w, 1 - sum(w)))
    sum <- 0
    for (i in 1:N) {
      # sum <- sum + log(exp(lpd_point[i, ]) %*% w_full)
      sum <- sum + dsmovetools::log_sum(lpd_point[i,] + w_full)
    }
    message(-as.numeric(sum))
    return(-as.numeric(sum))
  }
  gradient <- function(w) {
    stopifnot(length(w) == K - 1)
    w_full <- c(w, 1 - sum(w))
    grad <- rep(0, K - 1)
    for (k in 1:(K - 1)) {
      for (i in 1:N) {
        # grad[k] <- grad[k] + (exp_lpd_point[i, k] - exp_lpd_point[i, K])/(exp_lpd_point[i, ] %*% w_full)
        if(lpd_point[i, k] > lpd_point[i, K]) {
          grad[k] <- grad[k] + exp(
            c_log_subtract(lpd_point[i, k], lpd_point[i, K]) - dsmovetools::log_sum(lpd_point[i, ] + w_full)
          )
        } else {
          grad[k] <- grad[k] - exp(
            c_log_subtract(lpd_point[i, K], lpd_point[i, k]) - dsmovetools::log_sum(lpd_point[i, ] + w_full)
          )
        }
        if(any(!is.finite(grad))) {
          browser()
        }
      }
    }
    return(-grad)
  }
  ui <- rbind(rep(-1, K - 1), diag(K - 1))
  ci <- c(-1, rep(0, K - 1))
  o <- constrOptim(theta = rep(1/K, K - 1), f = negative_log_score_loo, 
                   grad = gradient, ui = ui, ci = ci, method = optim_method, 
                   control = optim_control)
  if(o$convergence != 0) {
    browser()
  }
  message(o$value)
  w <- o$par
  wts <- structure(c(w, 1 - sum(w)), names = paste0("model", 
                                                    1:K), class = c("stacking_weights"))
  return(wts)
}

wts = stacking_weights_local(
  lpd_point = lpd_point, 
  optim_control = list(reltol=1e-10, maxit = 1e3), 
  optim_method = 'BFGS'
)

# store the weights alongside their model replicate number
wts_df = data.frame(
  rep = str_extract(string = f, pattern = '[0-9]+'),
  w = as.numeric(wts)
)

ofile = file.path('output', 'mcmc', 'chain_weights', 'stacking_weights.rds')
saveRDS(wts_df, file = ofile)
