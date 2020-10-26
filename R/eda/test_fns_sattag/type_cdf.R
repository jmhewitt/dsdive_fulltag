type_cdf = function(pre_post_list) {

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
  
  summary_stats = df %>% 
    # merge dive classifications with observations
    dplyr::left_join(diveTypes, by = 'diveIds') %>% 
    # process pre/post observations separately
    dplyr::group_by(type) %>% 
    # sequential proportion of time in shallow periods
    dplyr::mutate(cum_shallow_obs = cumsum(diveType == 'Shallow') / n()) %>% 
    # return to normal
    ungroup()
  
  # reshape pre/post summary statistics
  summary_stats_wider = data.frame(
    pre = summary_stats %>% dplyr::filter(type == 'pre') %>% 
      dplyr::select(cum_shallow_obs) %>% unlist(),
    post = summary_stats %>% dplyr::filter(type == 'post') %>% 
      dplyr::select(cum_shallow_obs) %>% unlist()
  )
  
  # Cramer-von Mises like test: SSE for differences in summary statistics
  sum(apply(summary_stats_wider, 1, diff)^2)
}
