raw_data_plan = drake_plan(
  
  # seconds between depth observations
  sattag_timestep = 300,
  
  # CEE start times
  cee_starts = mdy_hms(c('5/15/19 15:04:50', '6/7/19 15:13:09', 
                         '8/6/19 17:37:43', '8/19/19 19:11:00'), tz = 'UTC'),
  
  # location of sattag series data
  depth_files = file_in(!!dir(path = file.path('data', 'raw'), 
                              pattern = 'series_', full.names = TRUE)),
  
  # location of message files associated with sattag series data
  message_files = file_in(!!dir(path = file.path('data', 'raw'), 
                                pattern = 'seriesrange_', full.names = TRUE)),
  
  # sex information for tags
  tag_sex = read.csv(file_in(!!file.path('data', 'raw', 'tag_sex.csv')), 
                     colClasses = 'factor'),
  
  # threshold for deep dives
  deep_dive_depth = 800
  
)