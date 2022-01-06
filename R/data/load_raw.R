load_raw = function(depth_files, template_bins, tag_info, dive_labels,
                    deep_depth_threshold, lon_lat_mean, timestep) {
  
  mapply(function(depth_file, labels) {
    
    # load data
    d = read.csv(file.path(depth_file))
    d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', tz = 'UTC')
    
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
    
    # get tag name
    tag_id = str_extract(depth_file, 'ZcTag[0-9]+')
    
    # extract cee start time
    exposure_time = (tag_info %>% 
      dplyr::filter(deployid == tag_id) %>% 
      dplyr::select(cee_start))[1,]
    
    # extract time at which baseline observations end
    baseline_end = ifelse(
      is.na(exposure_time), 
      tag_end + 1,
      (tag_info %>% 
         dplyr::filter(deployid == tag_id) %>% 
         dplyr::select(baseline_end))[1,]
    )
    
    # label observations as pre or post exposure
    exposed = as.numeric(d$Date >= exposure_time)
    exposed[is.na(exposed)] = 0
    baseline = as.numeric(d$Date < baseline_end)
    
    # add dive segmentation to depths time series
    d$diveId = labels$labels
    
    # classify dives by max observed depth
    diveTypes = d %>%
      dplyr::group_by(diveId) %>%
      dplyr::summarise(
        diveType = ifelse(any(diveId == 0), 'Unknown', ifelse(
          max(depth.standardized) > deep_depth_threshold, 'Deep', 'Shallow'
        ))
      )
    
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
      # dive segmentation
      diveIds = d$diveId,
      diveTypes = diveTypes,
      # exposure information
      exposure_time = exposure_time,
      exposed = exposed,
      baseline = baseline,
      baseline_end = baseline_end,
      # proximal location, near which data are collected
      proximal_loc = lon_lat_mean
    )
  }, depth_files, dive_labels, SIMPLIFY = FALSE)
  
}