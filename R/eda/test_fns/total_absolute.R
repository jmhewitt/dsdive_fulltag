total_absolute = function(pre_post_list) {
  # SSE for
  sum(diff(
    sapply(pre_post_list, function(d) {
      # total upward distance
      dx = diff(d$depths)
      sum(abs(dx))
    })
  )^2)
}
