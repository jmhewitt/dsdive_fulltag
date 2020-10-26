decreasing_type_seq = function(pre_post_list) {

  # merge data, for easier processing
  df = rbind(
    cbind(pre_post_list$pre, type = 'pre'),
    cbind(pre_post_list$post, type = 'post')
  ) 
  
  # classify dives as deep/shallow based on max observed depth
  diveTypes = df %>% 
    dplyr::group_by(diveIds) %>% 
    dplyr::summarise(diveType = factor(
      ifelse(max(depths) > 800, 'Deep', 'Shallow')
    ))
  
  # compute pre/post summary statistics
  summary_stats = df %>% 
    # merge dive classifications with observations
    dplyr::left_join(diveTypes, by = 'diveIds') %>% 
    # process pre/post observations separately
    dplyr::group_by(type) %>%  
    # Wilcoxon-style statistic to study sequence of dive types
    dplyr::summarise(W = sum(1:n() * as.numeric(diveType)))
  
  # SSE for differences in summary statistics
  diff(rev(summary_stats$W))
}
