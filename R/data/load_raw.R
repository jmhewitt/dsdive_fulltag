load_raw = function(depth_files, template_bins, tag_timelines,
                    deep_depth_threshold, lon_lat_mean, timestep) {
  
  mapply(function(depth_file) {
    
    # get tag name
    tag_id = str_extract(depth_file, 'ZcTag[0-9]+')
    
    # load data
    d = read.csv(file.path(depth_file))
    d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', tz = 'UTC')
    
    # validate all series data have unique timestamps
    d$dtime = c(0, diff(d$Date))
    if(any(diff(d$Date) == 0)) {
      if(tag_id == 'ZcTag111') {
        # remove duplicate data records from gonio source
        valid_source_rows = !((d$dtime == 0) & (d$original == 'gonio'))
      } else if(tag_id == 'ZcTag126') {
        valid_source_rows = 1:nrow(d)
        # remove duplicate data records from portal source
        valid_source_rows[
          which(((d$dtime == 0) & (d$original == 'gonio'))) - 1
        ] = FALSE
      } else {
        stop(paste('Duplicate series timestamps for ', tag_id))
      }
      d = d[valid_source_rows, ]
    } else {
      valid_source_rows = 1:nrow(d)
    }
    
    # forward difference: time elapsed from one observation to the next
    time_increments = difftime(
      time2 = d$Date[1:(length(d$Date)-1)], 
      time1 = d$Date[2:length(d$Date)], 
      units = 'secs'
    )
    
    # identify (inclusive) time ranges for data gaps
    d$gap_after = c(time_increments > timestep, FALSE)

    # map all depths to standardized bins
    d$depth.bin = sapply(d$Depth, function(depth) {
      which.min(abs(depth - template_bins$center))
    })
    d$depth.standardized = template_bins$center[d$depth.bin]
    
    # determine tag's temporal support
    tag_start = d$Date[1]
    tag_end = d$Date[length(d$Date)]
    
    # extract tag's timeline
    timeline = tag_timelines[[pmatch(tag_id, names(tag_timelines))]]
    
    # extract cee start times
    exposure_times = timeline$cee_segments$start
    
    # extract times at which baseline observations end
    baseline_ends = timeline$baseline_segments$end
    
    # label exposed observations
    if(nrow(timeline$cee_segments) > 0) {
      exposed = rowSums(
        do.call(cbind, lapply(1:nrow(timeline$cee_segments), function(i) {
          (timeline$cee_segments$start[i] <= d$Date) & 
            (d$Date <= timeline$cee_segments$end[i])
        }))
      ) > 0 
    } else {
      exposed = rep(FALSE, length(d$Date))
    }
    
    # label baseline observations
    baseline = rowSums(
      do.call(cbind, lapply(1:nrow(timeline$baseline_segments), function(i) {
        (timeline$baseline_segments$start[i] <= d$Date) & 
        (d$Date <= timeline$baseline_segments$end[i])
      }))
    ) > 0 
    
    # # add dive segmentation to depths time series
    # d$diveId = labels$labels[valid_source_rows]
    # 
    # # classify dives by max observed depth
    # diveTypes = d %>%
    #   dplyr::group_by(diveId) %>%
    #   dplyr::summarise(
    #     diveType = ifelse(any(diveId == 0), 'Unknown', ifelse(
    #       max(depth.standardized) > deep_depth_threshold, 'Deep', 'Shallow'
    #     ))
    #   )
    
    # enrich data with celestial covariates
    d$daytime = daytime(date = d$Date, lat = lon_lat_mean['lat'], 
                        lon = lon_lat_mean['lon'])
    d$moonlit = moonlit(date = d$Date, lat = lon_lat_mean['lat'], 
                        lon = lon_lat_mean['lon'])
    
    # package results
    list(
      # get name of tagged individual
      tag = tag_id,
      # data gap markers
      gap_after = d$gap_after,
      # standardized depths and times
      depth.bin = d$depth.bin,
      depths = d$depth.standardized,
      times = d$Date,
      # celestial covariates
      daytime = d$daytime,
      moonlit = d$moonlit,
      # # dive segmentation
      # diveIds = d$diveId,
      # diveTypes = diveTypes,
      # exposure information
      exposed = exposed,
      baseline = baseline,
      timeline = timeline,
      # proximal location, near which data are collected
      proximal_loc = lon_lat_mean
    )
  }, depth_files, SIMPLIFY = FALSE)
  
}