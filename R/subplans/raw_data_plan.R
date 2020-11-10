raw_data_plan = drake_plan(
  
  # seconds between depth observations
  sattag_timestep = 300,
  
  # location of dtag data
  dtag_files = file_in(!!dir(path = file.path('data', 'raw'),
                             pattern = 'aprh', full.names = TRUE)),
  
  # location of sattag series data
  depth_files = file_in(!!dir(path = file.path('data', 'raw'),
                              pattern = 'series_', full.names = TRUE)),
  
  # location of message files associated with sattag series data
  message_files = file_in(!!dir(path = file.path('data', 'raw'), 
                                pattern = 'seriesrange_', full.names = TRUE)),
  
  # sex and CEE information for tags
  tag_info = read.csv(file_in(!!file.path('data', 'raw', 'tag_info.csv')), 
                      colClasses = 'factor') %>% 
    dplyr::mutate(
      baseline_end = parse_date_time(x = baseline_end, orders = 'mdy HMS'),
      cee_start = parse_date_time(x = cee_start, orders = 'mdy HMS')
    ),

  # threshold for deep dives
  deep_dive_depth = 800
  
)