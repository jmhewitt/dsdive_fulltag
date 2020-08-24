# create a vector to label dive stages
stagevec = function (length.out, breaks) {
  rep(1:3, c(breaks[1] - 1, breaks[2] - breaks[1], length.out + 
               1 - breaks[2]))
}

# use 85% max depth rule to determine time in stages
times.stages = function(dives.obs) {
  do.call(rbind, lapply(dives.obs, function(d) {
    # extract depths
    depths = d$depth.bins$center[d$dive$depths]
    # find max depth, and stage threshold
    max.depth = max(depths)
    stage.thresh = .85 * max.depth
    # compute observed stage vector
    bottom.range = range(which(depths >= stage.thresh))
    if(length(unique(bottom.range))==1) {
      bottom.range[2] = bottom.range[1] + 1
    }
    stages = stagevec(length.out = length(depths), breaks = bottom.range)
    # linearly interpolate to find stage transition times
    t.inds = which(diff(stages)==1)
    t.stages = sapply(t.inds, function(ind) {
      # get start and end times/depths
      d0 = depths[ind]
      t0 = d$dive$times[ind]
      df = depths[ind+1]
      tf = d$dive$times[ind+1]
      # compute time at which stage.thresh is crossed
      if(df==d0) {
        mean(c(t0,tf))
      } else {
        (stage.thresh - d0)/(df-d0) * (tf-t0) + t0
      }
    })
    # return results
    data.frame(sub.time.min = t.stages[1]/60, 
               bottom.time.min = diff(t.stages)/60)
  }))
}