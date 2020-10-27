## Load your packages, e.g. library(drake).
source("./packages.R")

## Load your R files and subplans
invisible(lapply(list.files("./R", full.names = TRUE, recursive = TRUE), source))

source('R/plan.R')

# make(the_plan, lock_envir = FALSE, 
#      targets = grep('eda_tagsummary', the_plan$target, value = TRUE)[2])
options(clustermq.scheduler = "multicore")
make(the_plan, lock_envir = FALSE, 
     targets = grep('sattag_eda_plots', the_plan$target, value = TRUE),
     parallelism = 'clustermq', jobs = 5)


grep('sattag_eda_stats_12_1_conditional', the_plan$target, value = TRUE)


# make(the_plan, lock_envir = FALSE, targets = 'mcmc_samples_individual_nim_pkg_0')
# make(the_plan, lock_envir = FALSE, targets = 'posterior_report_single')


# make(the_plan, lock_envir = FALSE, 
#      targets = grep(pattern = 'exposure_eda', x = the_plan$target, value = TRUE))

# build graph components
graph = vis_drake_graph(the_plan, targets_only = TRUE)

# view graph
visNetwork::visHierarchicalLayout(graph, direction = "LR",
                                  edgeMinimization = FALSE)
