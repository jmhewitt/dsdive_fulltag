data_targets = list(
  
  # nominal location where data is collected
  tar_target(cape_hatteras_loc, c('lon' = -75.54, 'lat' = 35.23)),
  
  # seconds between depth observations
  tar_target(sattag_timestep, 300),
  
  # deep dive threshold (m)
  tar_target(deep_dive_depth, 800),
  
  
  tar_target(sattag_dir, file.path('data', 'sattag')),
  
  # sattag series data
  tar_target(
    name = depth_files, 
    command = dir(path = sattag_dir, pattern = 'Series\\.csv', 
                  full.names = TRUE),
    format = 'file'
  ),
  
  # sattag messages
  tar_target(
    name = message_files, 
    command = dir(
      path = sattag_dir, pattern = 'SeriesRange\\.csv', full.names = TRUE
    ),
    format = 'file'
  ),
  
  # formats of which tag metadata are stored in
  tar_target(
    name = date_formats, 
    command = c('mdy IMS p', 'ymd HMS')
  ),
  
  # tag metadata: sex and CEE info
  tar_target(
    name = tag_info, 
    command = read.csv(file.path(sattag_dir, 'tag_info.csv'), 
                       colClasses = 'factor') %>%
        dplyr::mutate(
          baseline_end = parse_date_time(
            x = baseline_end, orders = date_formats, tz = 'UTC'
          ),
          cee_start = parse_date_time(
            x = cee_start, orders = date_formats, tz = 'UTC'
          )
        )
  ),
  
  # load tags into standard format
  tar_target(
    name = raw_sattags,
    command = load_raw(depth_files = depth_files, 
                       template_bins = template_bins, tag_info = tag_info, 
                       dive_labels = dive_labels, deep_depth_threshold = 800, 
                       lon_lat_mean = cape_hatteras_loc,
                       timestep = sattag_timestep)
  )
  
)
