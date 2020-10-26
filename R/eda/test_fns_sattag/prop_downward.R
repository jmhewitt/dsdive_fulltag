prop_downward = function(pre_post_list) {
  # SSE for
  sum(diff(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # proportion of downward transitions
      mean( sign(diff(d$depth_bins)) == 1 )
    })
  )^2)
}
