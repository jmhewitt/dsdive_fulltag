decreasing_total_downward = function(pre_post_list) {
  diff(rev(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # total upward distance
      dx = diff(d$depths)
      sum(dx[dx<0])
    })
  ))
}
