load_raw = function(depth_files, template_bins, tag_info, dive_labels) {
  
  lapply(1:length(depth_files), function(i) {
    
    # load data
    d = read.csv(file.path(depth_files[i]))
    d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', tz = 'UTC')
    
    # map all depths to standardized bins
    d$depth.bin = sapply(d$Depth, function(depth) {
      which.min(abs(depth - template_bins$center))
    })
    d$depth.standardized = template_bins$center[d$depth.bin]
    
    # determine tag's temporal support
    tag_start = d$Date[1]
    tag_end = d$Date[length(d$Date)]
    
    # get tag name
    tag_id = str_extract(depth_files[i], 'ZcTag[0-9]+')
    
    # determine which CEE's overlap with tag record
      cees_experienced = tag_info %>% 
        dplyr::filter(deployid == tag_id) %>% 
        dplyr::select(cee_start) %>%
      unlist()
      
    # pre/post exposure and baseline labels
    if(!is.na(cees_experienced) == 0) {
      exposure_time = Inf
      exposed = numeric(nrow(d))
      baseline = numeric(nrow(d)) + 1
      baseline_end = tag_end
    } else {
      exposure_time = as.POSIXct(
        cees_experienced, origin = '1970-01-01 00:00.00 UTC'
      )
      baseline_end = as.POSIXct(
        tag_info %>% 
          dplyr::filter(deployid == tag_id) %>% 
          dplyr::select(baseline_end) %>% 
          unlist(), 
        origin = '1970-01-01 00:00.00 UTC'
      )
      exposed = as.numeric(d$Date >= exposure_time)
      baseline = as.numeric(d$Date < baseline_end)
    }
    
    # add dive segmentation to depths time series
    labelInd = which(
      sapply(dive_labels, function(x) x$tag) == as.character(d$DeployID[1])
    )
    d$diveId = dive_labels[[labelInd]]$labels

    # classify dives by max observed depth
    diveTypes = d %>%
      dplyr::group_by(diveId) %>%
      dplyr::summarise(diveType = ifelse(max(depth.standardized) > 800,
                                         'Deep', 'Shallow'))
    
    # package results
    list(
      # get name of tagged individual
      tag = levels(d$DeployID)[d$DeployID[1]],
      # standardized depths and times
      depth.bin = d$depth.bin,
      depths = d$depth.standardized,
      times = d$Date,
      # dive segmentation
      diveIds = d$diveId,
      diveTypes = diveTypes,
      # exposure information
      exposure_time = exposure_time,
      exposed = exposed,
      baseline = baseline,
      baseline_end = baseline_end
    )
  })
   
}