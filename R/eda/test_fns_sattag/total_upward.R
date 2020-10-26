total_upward = function(pre_post_list) {
  # SSE for
  sum(diff(
    sapply(pre_post_list[c('pre', 'post')], function(d) {
      # total upward distance
      dx = diff(d$depths)
      sum(dx[dx>0])
    })
  )^2)
}
