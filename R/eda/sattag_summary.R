sattag_summary = function(dat) {
  dat %>% dplyr::summarise(
    'prop_surface' = mean(depth_bins == 1),
    'prop_downward' = mean(descending),
    'prop_upward' = mean(ascending),
    'prop_no_change' = mean(no_change),
    'total_absolute' = sum(abs(ddepths)),
    'total_downward' = sum(ddepths * descending),
    'total_upward' = sum(ddepths * ascending),
    'depth_seq' = sum(1:n() * depths),
    'type_seq' = sum(1:n() * as.numeric(factor(diveType)))
  ) 
}