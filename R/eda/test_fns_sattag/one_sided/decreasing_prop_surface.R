decreasing_prop_surface = function(pre_post_list) {
  diff(rev(
    # proportion of time spent in shallowest depth bin
    sapply(
      pre_post_list[c('pre', 'post')], 
      function(d) mean(d$depth_bins == 1) 
    )
  ))
}
