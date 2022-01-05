library(targets)
library(future)
library(future.batchtools)

# basic check for existence of SLURM job submission command to determine if 
# running on a SLURM-enabled server
if(system('command -v sbatch') == 0) {
  plan(batchtools_slurm, template = file.path("hpc", "slurm_batchtools.tmpl"))
} else {
  plan(multisession)
}

# set packages to load
tar_option_set(
  packages = c('dplyr', 'lubridate', 'ggplot2', 'ggthemes', 'stringr', 
               'nimble', 'expm', 'pryr', 'suncalc', 'tarchetypes', 'coda',
               'tidyr', 'future', 'future.batchtools', 'viridis', 'splines2'),
  deployment = 'main'
)

## load R files and workflows
lapply(list.files("R", full.names = TRUE, recursive = TRUE, pattern = '\\.R'), 
       source)

# assemble workflow
c(
  dir_targets,
  data_targets, 
  depth_bin_imputation_targets,
  dive_segmentation_targets,
  model_discretization_target,
  nimble_targets,
  validation_targets,
  cee_targets,
  eda_targets,
  post_sim_targets,
  plot_targets
)