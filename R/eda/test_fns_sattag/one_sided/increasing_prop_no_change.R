increasing_prop_no_change = function(pre_post_list) {
  diff(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # proportion of transitions where no change was seen
      mean( sign(diff(d$depth_bins)) == 0 )
    })
  )
}
