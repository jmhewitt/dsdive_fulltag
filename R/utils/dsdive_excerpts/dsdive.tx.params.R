#' Compute transition parameters for dive trajectories across discrete depths
#'
#' Computes elements of a transition matrix when given model parameters and 
#' time/space locations for a dive model that has 3 stages, PRIMARY DESCENT 
#' (PD), INTERMEDIATE BEHAVIORS (IB), and PRIMARY ASCENT (PA).
#'   
#' @param t0 time at which transition parameters should be computed
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param d0 the depth bin at which transition parameters should be computed
#' @param d0.last the previous depth bin in which the trajectory was.  If 
#'   \code{NULL}, then the autoregressive component will be skipped.
#' @param s0 the dive stage (PRIMARY DESCENT==1, INTERMEDIATE BEHAVIORS==2, 
#'   PRIMARY ASCENT==3) for which transition parameters should be computed.  
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the PRIMARY DESCENT 
#'  (PD), INTERMEDIATE BEHAVIORS (IB), and PRIMARY ASCENT (PA) dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the PD, IB, and PA stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the IB stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the PA stage at the next depth transition
#' @param t0.dive Time at which dive started
#' @param shift.params If not \code{NULL}, then will apply a directional bias 
#'   to depth bin transition probabilities.  \code{shift.params} should be a
#'   vector \code{c(df, p.bias, tf, lambda.bias)}. The first element is the 
#'   target depth bin to travel to.
#'   The second element controls the strength of the bias.  The second element 
#'   should be a positive value, with 1.125 being the recommended value.  Using 
#'   \code{shift.params[2] = 1.125} will reinforce natural movement in good 
#'   directions of travel, and replace "poor" directions of travel with bias in 
#'   "good" directions of travel.  Setting \code{shift.params[2] = 1} will 
#'   replace "poor" directions of travel with random walk proposals.  Values 
#'   less than 1 will push poor directions of travel closer toward random walks
#'   from their natural tendencies.
#' @param t.stage2 time at which second stage was entered.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#'   
#' @example examples/txparams.R
#' 
#' @importFrom stats plogis
#' 
#' @export
#' 
dsdive.tx.params = function(depth.bins, d0, s0, beta, lambda) {
  
  num.depths = nrow(depth.bins)
  
  # extract transition rate
  rate = lambda[s0] / (2 * depth.bins[d0, 2])
  
  # define transition probabilities between dive bins
  if(d0 == 1) {
    
    # surface can only transition downward
    nbrs = 2
    
    # can only leave surface depth bin in stages 1 or 2
    probs = ifelse(s0 < 3, 1, 0)
    
  } else if(d0 == 2) {
    
    nbrs = c(1,3)
    
    # can only return to surface bin in stage 3
    probs.down = ifelse(s0==3, beta[2], 1)
    probs = c(1-probs.down, probs.down)
    
  } else if(d0 == num.depths) {
    
    # bottom bin is a boundary; transition is deterministic
    nbrs = num.depths - 1
    probs = 1
    
  } else {
    
    nbrs = d0 + c(-1, 1)
    
    probs.down = ifelse(s0==1, beta[1], ifelse(s0==2, .5, beta[2]))
    probs = c(1-probs.down, probs.down)
    
  }
  
  # package results
  list(rate = rate, probs = probs, labels = nbrs)
}