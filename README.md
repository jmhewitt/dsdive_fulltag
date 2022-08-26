- Run make.R to prepare data for model fitting
- Run Rscripts/model_fitting/dothescience.job on HPC to fit model with 
  multiple random starts, and a single fixed start
- Run additional scripts to generate additional posterior output
- Run additional scripts to generate figures and tables

The scripts are not all done in a single workflow because parallelism is useful,
and we deal with this differently.  Note that we can change the script make.R
so that it automatically kicks off the additional scripts to save some time, 
but, of course, we won't have the final scripts be able to be automatically run.
Maybe, though, this is fine b/c we want to work on automating the generation of 
main model output, and we can deal with all of the additional post-processing 
separately.