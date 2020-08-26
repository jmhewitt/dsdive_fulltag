## Load your packages, e.g. library(drake).
source("./packages.R")

## Load your R files and subplans
invisible(lapply(list.files("./R", full.names = TRUE, recursive = TRUE), source))

source('R/plan.R')

make(the_plan, lock_envir = FALSE)

# build graph components
graph = vis_drake_graph(the_plan, targets_only = TRUE)

# view graph
visNetwork::visHierarchicalLayout(graph, direction = "LR",
                                  edgeMinimization = FALSE)

r_make()
