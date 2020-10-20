mean_upward_velocity = function(pre_post_list) {
  # SSE for
  sum(diff(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # average upward velocity
      dpdt = diff(d$p) / diff(as.numeric(d$times))
      mean(dpdt[dpdt > 0])
    })
  )^2)
}
