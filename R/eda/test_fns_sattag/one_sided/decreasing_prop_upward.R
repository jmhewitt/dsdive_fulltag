decreasing_prop_upward = function(pre_post_list) {
  diff(rev(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # proportion of upward transitions
      mean( sign(diff(d$depth_bins)) == -1 )
    })
  ))
}
