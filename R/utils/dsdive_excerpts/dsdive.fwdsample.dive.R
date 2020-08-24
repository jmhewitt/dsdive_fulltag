#' Simulate dive trajectories across discrete depth bins
#'
#' The method will simulate dive trajectories from initial conditions until the 
#' trajectory is observable at \code{tf}, or a maximum number of transitions 
#' has been exceeded.  The dive simulation is bridged, so the trajectory will
#' also stop diving after returning to the surface.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param d0 the depth bin at which transition parameters should be computed
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the DIVING, SUBMERGED, 
#'  and SURFACING dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the DIVING, SUBMERGED, and SURFACING stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the SUBMERGED stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the SURFACING stage at the next depth transition
#' @param t0 time at which transition parameters should be computed
#' @param tf time at which sampling should end after
#' @param steps.max maximum number of transitions to sample before stopping, 
#'   regardless of whether \code{tf} is reached.
#' @param dur0 time spent at location \code{d0}.  If \code{NULL}, then a 
#'   duration in state \code{d0} will be sampled, otherwise a new state will 
#'   be sampled first, then sampling will continue from the new state at time 
#'   \code{t0 + dur0}.
#' @param nsteps if \code{nsteps} is not \code{NULL}, then the sampler will
#'   be reconfigured to sample exactly \code{nsteps} transitions instead of 
#'   attempting to sample until the trajectory is observable at time \code{tf}.
#' @param s0 dive stage in which forward simulation begins
#' @param t0.dive Time at which dive started
#' @param shift.params Optional arguments to bias sampling toward a specific 
#'   node.  See documentation for \code{dsdive.tx.params} for more detail.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' @param T1 stage 1 length (in same units as t0, d0, tf)
#' @param T2 stage 2 length (in same units as t0, d0, tf)
#' 
#' @return A \code{dsdive} object, which is a \code{list} with the following 
#'   vectors:
#'   \describe{
#'     \item{depths}{Record of which depth bin indices the trajectory visited}
#'     \item{durations}{Record of amount of time spent in each depth bin}
#'     \item{times}{The time at which each depth bin was entered}
#'     \item{stages}{The stage at which each depth bin was entered}
#'   }
#'   
#' 
#' @example examples/dsdive.fwdsample.dive.R
#' 
#' 
#' @export
#' 
#'
dsdive.fwdsample.dive = function(depth.bins, beta, lambda, t0, steps.max, 
                                 T1, T2) {
  
  #
  # sample stage 1 portion of dive
  #
  
  tf1 = t0 + T1
  d1 = dsdive.fwdsample.fixedstage(depth.bins = depth.bins, d0 = 1, 
                                   beta = beta, lambda = lambda, t0 = t0, 
                                   tf = tf1, steps.max = steps.max, s0 = 1)
  
  # adjust dive s.t. we stop stage 1 dynamics when stage 1 ends
  n1 = length(d1$depths)
  d1$durations[n1] = tf1 - d1$times[n1]
  
  #
  # sample stage 2 portion of dive
  #
  
  tf2 = tf1 + T2
  d2 = dsdive.fwdsample.fixedstage(depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   d0 = d1$depths[n1],  
                                   dur0 = NULL,
                                   t0 = tf1, 
                                   tf = tf2, steps.max = steps.max, s0 = 2)
  
  
  # adjust dive s.t. we stop stage 2 dynamics when stage 2 ends
  n2 = length(d2$depths)
  d2$durations[n2] = tf2 - d2$times[n2]
  
  #
  # sample stage 3 portion of dive
  #
  
  d3 = dsdive.fwdsample.fixedstage(depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   d0 = d2$depths[n2],  
                                   dur0 = NULL,
                                   t0 = tf2, 
                                   tf = Inf, steps.max = steps.max, s0 = 3)
  
  #
  # package dive
  #
  
  res = list(
    depths = c(d1$depths, d2$depths, d3$depths),
    times = c(d1$times, d2$times, d3$times),
    stages = c(d1$stages, d2$stages, d3$stages)
  )
  
  res$durations = c(diff(res$times), Inf)
  
  class(res) = 'dsdive'
  
  res
}