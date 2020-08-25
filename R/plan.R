lapply(list.files("./R/subplans", full.names = TRUE, recursive = TRUE), source)

the_plan = bind_plans(
  raw_data_plan, 
  depth_bin_imputation_plan,
  dive_segmentation_plan,
  nimble_plan
  #report_plan,
  #validation_plan,
  #simulation_plan
)