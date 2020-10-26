increasing_total_absolute = function(pre_post_list) {
  diff(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # total upward distance
      dx = diff(d$depths)
      sum(abs(dx))
    })
  )
}
