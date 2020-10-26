decreasing_total_absolute = function(pre_post_list) {
  diff(rev(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # total upward distance
      dx = diff(d$depths)
      sum(abs(dx))
    })
  ))
}
