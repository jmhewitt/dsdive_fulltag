library(targets)

# set packages to load
tar_option_set(
  packages = c('dplyr', 'lubridate', 'ggplot2', 'ggthemes', 'stringr', 
               'nimble', 'expm', 'pryr', 'suncalc')
)

## load R files and workflows
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

# assemble workflow
c(
  dir_targets,
  data_targets, 
  depth_bin_imputation_targets,
  dive_segmentation_targets,
  model_discretization_target,
  nimble_targets
)