prop_surface = function(pre_post_list) {
  # SSE for
  sum(diff(
    # proportion of time spent in shallowest depth bin
    sapply(
      pre_post_list[c('pre', 'post')], 
      function(d) mean(d$depth_bins == 1) 
    )
  )^2)
}
