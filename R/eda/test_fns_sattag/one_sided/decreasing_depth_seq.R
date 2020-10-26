decreasing_depth_seq = function(pre_post_list) {

  # merge data, for easier processing
  df = rbind(
    cbind(pre_post_list$pre, type = 'pre'),
    cbind(pre_post_list$post, type = 'post')
  ) 
  
  # compute pre/post summary statistics
  summary_stats = df %>% 
    # process pre/post observations separately
    dplyr::group_by(type) %>%  
    # Wilcoxon-style statistic to study sequence of dive types
    dplyr::summarise(W = sum(1:n() * depths))
  
  # SSE for differences in summary statistics
  diff(rev(summary_stats$W))
}
