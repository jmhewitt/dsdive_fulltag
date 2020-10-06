bin_counts = function(pre_post_list) {

  df = rbind(
    cbind(pre_post_list$pre, type = 'pre'),
    cbind(pre_post_list$post, type = 'post')
  )
  
  # SSE for differences in marginal depth bin distributions
  sum(apply(table(df %>% dplyr::select(depths, type)), 1, diff)^2)
}
