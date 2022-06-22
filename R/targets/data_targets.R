data_targets = list(
  
  # nominal location where data is collected
  tar_target(cape_hatteras_loc, c('lon' = -75.54, 'lat' = 35.23)),
  
  # seconds between depth observations
  tar_target(sattag_timestep, 300),
  
  # deep dive threshold (m)
  tar_target(deep_dive_depth, 800),
  
  # series tag directory
  tar_target(
    name = sattag_dir, 
    command = file.path('data', 'scaled_source_ms-dev-main', 'data', 
                        '03_SERIES_SERIESRANGE_ALL', 'L2_gonio_pressuresensor')
  ),
  
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
  
  # use CEE metadata to identify baseline and exposed periods in data
  tar_target(
    name = tag_timelines,
    command = {
      
      # load cee metadata, parse times, define return-to-baseline via 24h rule
      cee_info = read.csv(
        file.path('data', 'scaled_source_ms-dev-main', 'data', 
                  'tagstreams_by_cee_details_filt.csv')
      ) %>% mutate(
        tstart = as.POSIXct(strptime(
          x = paste(date_YYYYMMDD, start_HHMMSS),
          format = '%Y%m%d %H%M%S', 
          tz = 'UTC'
        )),
        tend = as.POSIXct(strptime(
          x = paste(date_YYYYMMDD, end_HHMMSS),
          format = '%Y%m%d %H%M%S', 
          tz = 'UTC'
        )),
        treturn_to_baseline = tend + 24 * 3600
      )
      
      # enumerate baseline and cee periods for each tag
      timelines = lapply(depth_files, function(f) {
        
        # load series data, parse times
        seriesdata = read.csv(f) %>% 
          mutate(
            Date = as.POSIXct(Date, origin = '1970-01-01 00:00.00 UTC', 
                              tz = 'UTC')
          )
        
        # extract tag name
        tag = unique(seriesdata$DeployID)
        
        # basic data validation
        if(length(tag) > 1) {
          stop('Series data file contains measurements from multiple tags.')
        }
        
        # tag deployment period
        tag_range = range(seriesdata$Date)
        
        # cee's used to split data into distinct periods
        cees = cee_info %>% 
          filter(
            # only focus on cee's associated with current tag
            deployid == tag,
            # CEEs must begin during data collection
            tag_range[1] <= tstart, 
            tstart <= tag_range[2]
          )
        
        # build timeline of baseline/exposed periods for tag
        critical_times = c(tag_range, cees$tstart, cees$treturn_to_baseline)
        baseline_end = c(0,2, rep(1, nrow(cees)), rep(0, nrow(cees)))
        
        # remove timeline points past the tag range, such as a return to 
        # baseline after tag ends
        observable_times = (tag_range[1] <= critical_times) & 
          (critical_times <= tag_range[2])
        critical_times = critical_times[observable_times]
        baseline_end = baseline_end[observable_times]
        
        # sort timeline
        o = order(critical_times)
        critical_times = critical_times[o]
        baseline_end = baseline_end[o]
        
        # validate timeline
        if(length(unique(critical_times)) != length(critical_times)) {
          stop(paste('Duplicate times in timeline for', tag))
        }
        
        #
        # enumerate baseline periods
        #
        
        # baseline period starts
        start_inds = which(baseline_end == 0) 
        # only keep baseline periods that have clear ends
        keep = baseline_end[start_inds + 1] %in% 1:2
        start_inds = start_inds[keep]
        
        # basic data validation
        if(!all(keep)) {
          stop(paste('Potentially inconsistent baseline timeline for', tag))
        }
        
        # assemble baseline periods
        baseline_segments = data.frame(
          start = critical_times[start_inds],
          end = critical_times[start_inds + 1]
        ) %>% mutate(
          duration_days = (as.numeric(end) - as.numeric(start)) / 3600 / 24
        )
        
        #
        # enumerate cee periods
        #
        
        # cee period starts
        start_inds = which(baseline_end == 1)
        # only keep cee periods that have clear ends
        keep = baseline_end[start_inds + 1] %in% c(0,2)
        start_inds = start_inds[keep]
        
        # basic data validation
        if(!all(keep)) {
          stop(paste('Potentially inconsistent CEE timeline for', tag))
        }
        
        # assemble cee periods
        cee_segments = data.frame(
          start = critical_times[start_inds],
          end = critical_times[start_inds + 1]
        ) %>% mutate(
          duration_days = (as.numeric(end) - as.numeric(start)) / 3600 / 24
        )
        
        # package results
        list(
          tag = tag,
          baseline_segments = baseline_segments,
          cee_segments = cee_segments
        )
      })
      
      # label and return output
      names(timelines) = sapply(timelines, function(x) x$tag)
      timelines
    }
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
