## Load your packages, e.g. library(drake).
source("./packages.R")

## Load your R files and subplans
invisible(lapply(list.files("./R", full.names = TRUE, recursive = TRUE), source))

source('R/plan.R')

make(the_plan, lock_envir = FALSE)