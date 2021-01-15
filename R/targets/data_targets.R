data_targets = list(
  
  # seconds between depth observations
  tar_target(sattag_timestep, 300),
  
  # deep dive threshold (m)
  tar_target(deep_dive_depth, 800),
  
  # sattag series data
  tar_target(
    name = depth_files, 
    command = dir(path = sattag_dir, pattern = 'series_', full.names = TRUE),
    format = 'file'
  ),
  
  # sattag messages
  tar_target(
    name = message_files, 
    command = dir(path = sattag_dir, pattern = 'seriesrange_', full.names = TRUE),
    format = 'file'
  ),
  
  # tag metadata: sex and CEE info
  tar_target(
    name = tag_info, 
    command = read.csv(file.path(sattag_dir, 'tag_info.csv'), 
                       colClasses = 'factor') %>%
        dplyr::mutate(
          baseline_end = parse_date_time(
            x = baseline_end, orders = 'mdy HMS', tz = 'UTC'
          ),
          cee_start = parse_date_time(
            x = cee_start, orders = 'mdy HMS', tz = 'UTC'
          )
        )
  )
  
)
