#' Simulate observing a fully observed dive trajectory at specified timepoints
#' 
#' Given the components of a completely observed dive trajectory, extract the 
#' depth and stage of a trajectory at the target times \code{t.obs}.
#' 
#' @param depths Complete record of depth bin indices visited
#' @param times Times at which each of \code{depths} was visited
#' @param stages Stages at which each of the \code{depths} was visited
#' @param t.obs Times at which the trajectory should be observed
#' 
#' @example examples/observe.R
#' 
#' @export
#' 
dsdive.observe = function(depths, times, stages, t.obs) {
  tind = findInterval(t.obs, times)
  
  res = list(
    depths = depths[tind],
    stages = stages[tind],
    times = t.obs
  )
  
  class(res) = 'dsobs'
  
  res
}