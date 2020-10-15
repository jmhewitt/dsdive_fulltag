load_raw_dtag = function(dtag_files, cee_starts) {
  
  lapply(1:length(dtag_files), function(i) {
    
    # extract tag name from filename
    tag = gsub(pattern = '(\\..*|aprh.*)', replacement = '', 
               x = basename(dtag_files[i]))
    
    # load metadata
    meta = xmlParse(
      file = dir(path = dirname(dtag_files[i]), 
                 pattern = paste(tag, '.*\\.xml', sep = ''), 
                 full.names = TRUE)
    )
    
    # extract tag start time from metadata
    tag_start = strptime(
      x = xpathSApply(doc = meta, path = '//EVENT[START/@STATE = "NEW"]/@TIME'), 
      format = '%Y,%m,%d,%H,%M,%S', tz = 'UTC'
    )
    
    # load tag data
    d = readMat(file.path(dtag_files[i]))
    
    # associate times with datapoints
    times = seq(from = tag_start, length.out = length(d$p), 
                by = as.numeric(1/d$fs))
    
    # determine tag's temporal support
    tag_end = tail(times, 1)
    
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
      exposed = times >= exposure_time
    }
    
    # add to results
    d$tag = tag
    d$exposed = exposed
    d$times = times
    d$exposure_time = exposure_time
    
    # return package
    d
  })
   
}