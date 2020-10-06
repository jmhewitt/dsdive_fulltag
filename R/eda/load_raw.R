load_raw = function(depth_files, template_bins, cee_starts) {
  
  lapply(depth_files, function(f) {
    
    # load data
    d = read.csv(file.path(f))
    d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', tz = 'UTC')
    
    # map all depths to standardized bins
    d$depth.bin = sapply(d$Depth, function(depth) {
      which.min(abs(depth - template_bins$center))
    })
    d$depth.standardized = template_bins$center[d$depth.bin]
    
    # determine tag's temporal support
    tag_start = d$Date[1]
    tag_end = d$Date[length(d$Date)]
    
    # determine which CEE's overlap with tag record
    cees_experienced = which(
      (tag_start <= cee_starts) & (cee_starts <= tag_end)
    )
    
    # pre/post exposure labels
    if(length(cees_experienced) == 0) {
      exposure_time = Inf
      exposed = numeric(length(d$depths))
    } else {
      exposure_time = cee_starts[min(cees_experienced)]
      exposed = as.numeric(d$Date >= exposure_time)
    }
    
    # package results
    list(
      # get name of tagged individual
      tag = levels(d$DeployID)[d$DeployID[1]],
      # standardized depths and times
      depth.bin = d$depth.bin,
      depths = d$depth.standardized,
      times = d$Date,
      # exposure information
      exposure_time = exposure_time,
      exposed = exposed
    )
  })
   
}