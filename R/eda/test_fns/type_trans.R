type_trans = function(pre_post_list) {

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
  
  # deep/shallow sequence in "pre" observations
  pre_seq = df %>% 
    dplyr::left_join(diveTypes, by = 'diveIds') %>% 
    dplyr::filter(type == 'pre') %>% 
    dplyr::select(diveType) %>% 
    unlist()
  
  # deep/shallow sequence in "post" observations
  post_seq = df %>% 
    dplyr::left_join(diveTypes, by = 'diveIds') %>% 
    dplyr::filter(type == 'post') %>% 
    dplyr::select(diveType) %>% 
    unlist()
  
  # SSE between contingency tables of Deep-shallow transitions
  sum((
    table(pre_seq[-length(pre_seq)], pre_seq[-1]) / length(pre_seq) - 
    table(post_seq[-length(post_seq)], post_seq[-1]) / length(pre_seq)
  )^2)
}
