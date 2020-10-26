decreasing_prop_downward = function(pre_post_list) {
  diff(rev(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # proportion of downward transitions
      mean( sign(diff(d$depth_bins)) == 1 )
    })
  ))
}
