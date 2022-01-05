library(targets)

tar_make_future(
  names = c('validate_deep_surv_distn'), 
  workers = 400
)

# 'cee_dive_response_summaries'

# TODO: 'deep_dive_time_preds_wrt_covariates' in post_sim_targets.R and 
# fwd_sim_to_depth_fixed_covs functions, but only if the validation results look
# fine.