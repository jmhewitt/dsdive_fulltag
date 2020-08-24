identify_dive_endpoints = function(depth_files, dive_labels, template_bins, 
                                   sattag_timestep) {
  
  endpoints = mapply(FUN = function(depth_file, labels) { 
  
    # extract tag name
    tag.name = labels$tag
    
    # extract dive labels
    labels = labels$labels
    
    # load data
    d = read.csv(file = depth_file)
    d$Date = anytime(paste(d$Day, d$Time))
    
    # map all depths to standardized bins
    d$depth.bin = sapply(d$Depth, function(depth) {
      which.min(abs(depth - template_bins$center))
    })
    d$depth.standardized = template_bins$center[d$depth.bin]
    
    #
    # initialize dive start/end random effect associations and priors
    #
    
    # determine conflict points that are adjacent to dives
    endpoint.inds = matrix(nrow = 0, ncol = 5)
    colnames(endpoint.inds) = c('ind', 't_lwr', 't_upr', 'dive_end', 
                                'dive_start')
    for(ind in which(labels==0)) {
      
      # verify conflict point is adjacent to at least one dive
      if(any(labels[ind + c(-1,1)] != 0)) {
        
        # get time associated with conflict point
        t_ind = as.numeric(d$Date[ind])
        
        # associate with end time for previous dive
        dive_end = ifelse(labels[ind - 1] != 0, labels[ind - 1], 0)
        
        # associate with start time for next dive
        dive_start = ifelse(labels[ind + 1] != 0, labels[ind + 1], 0)
        
        # specify prior distribution and dive associations for conflict point
        endpoint.inds = rbind(
          endpoint.inds, 
          c(ind, t_ind + sattag_timestep * c(-1,1), dive_end, dive_start)
        )
      }
      
    }
    
    #
    # merge random effects with overlapping priors
    #
    
    # convert to data.frame for easier manipulation
    endpoint.inds = data.frame(endpoint.inds)
    
    # merge random effects that have overlapping priors
    for(overlapping_ind in 
        which(endpoint.inds$t_lwr[2:nrow(endpoint.inds)] < 
              endpoint.inds$t_upr[1:(nrow(endpoint.inds)-1)])) {
      
      # merged row
      m = c(endpoint.inds$ind[overlapping_ind],
            endpoint.inds$t_lwr[overlapping_ind],
            endpoint.inds$t_upr[overlapping_ind + 1],
            endpoint.inds$dive_end[overlapping_ind],
            endpoint.inds$dive_start[overlapping_ind + 1])
      
      # duplicate row in list of priors
      endpoint.inds[overlapping_ind,] = m
      endpoint.inds[overlapping_ind + 1,] = m
      
    }
    
    # remove duplicate rows
    endpoint.inds = unique(endpoint.inds)
    
    # assign ids to endpoints
    endpoint.inds$endpoint.id = 1:nrow(endpoint.inds)
    
    # # verify each random effect is associated with at most one dive start/end
    # table(endpoint.inds$dive_end)[table(endpoint.inds$dive_end) > 1]
    # table(endpoint.inds$dive_start)[table(endpoint.inds$dive_start) > 1]
    
    #
    # associate individual observations with dives
    #
    
    # raw depth and observation time series
    obs = list(depths = d$depth.bin, times = as.numeric(d$Date))
    
    # indices of first potential observations for dives
    dive.ranges = data.frame(t(apply(endpoint.inds, 1, function(r) {
      start.ind = which(r['t_lwr'] == obs$times) + 1
      if(length(start.ind) == 0) {
        start.ind = NA
      }
      c(dive.id = as.numeric(r['dive_start']),
        start.ind = start.ind)
    }))) %>% dplyr::filter(dive.id > 0)
    
    # append indices of last potential observations for dives
    dive.ranges$end.ind = as.numeric(apply(
      dive.ranges %>% left_join(endpoint.inds, by = c('dive.id' = 'dive_end')), 
      1, 
      function(r) { which(r['t_upr'] == obs$times) - 1 }
    ))
    
    # remove dives that do not have fully-specified observation records
    dive.ranges = dive.ranges[complete.cases(dive.ranges),]
    
    # add indicator for whether dive is deep (type == 1) or shallow (type == 2)
    dive.ranges$type = factor(apply(dive.ranges, 1, function(r) {
      max_d = template_bins$center[
        max(d$depth.bin[r['start.ind']:r['end.ind']])
      ]
      ifelse(max_d > 800, 'Deep', 'Shallow')
    }))
    
    # add endpoint links to dive information
    dive.ranges = dive.ranges %>% 
      # link T0 endpoint
      left_join(endpoint.inds %>% mutate(T0_endpoint = endpoint.id) %>% 
                  dplyr::select(T0_endpoint, dive_start), 
                by = c(dive.id = 'dive_start')) %>% 
      # link T3 endpoint
      left_join(endpoint.inds %>% mutate(T3_endpoint = endpoint.id) %>% 
                  dplyr::select(T3_endpoint, dive_end),
                by = c(dive.id = 'dive_end'))
    
    #
    # (quality control) determine which dives to analyze
    #
    
    # extract all unique dive id's
    dive.ids = sort(unique(labels[labels>0]))
    
    # minimal flags to denote which dives will be analyzed
    dive.flags = sapply(dive.ids, function(id) {
      
      if(id %in% dive.ranges$dive.id) {
        
        # extract times associated with the dive
        rowInd = which(id == dive.ranges$dive.id)
        times = as.numeric(
          d$Date[seq(from = dive.ranges[rowInd, 'start.ind'],
                     to = dive.ranges[rowInd, 'end.ind'])]
        )
        
        ifelse(all(
          # dive has complete set of depth observations
          id %in% dive.ranges$dive.id, 
          all(diff(times) == sattag_timestep),
          # dive is associated with random effect for end-of-dive
          id %in% endpoint.inds$dive_end, 
          # dive is associated with random effect for start-of-dive
          id %in% endpoint.inds$dive_start), 
          TRUE, FALSE)
      } else {
        FALSE
      }
    })
    
    # collect processing data in a single location
    dives.processed = list(
      # define priors for dive start/end random effects, and associate with dives
      endpoint.inds = endpoint.inds,
      # associate observations and covariates to each dive
      dive.ranges = dive.ranges,
      # indicators for which dives to analyze
      dive.flags = dive.flags,
      # depth and time observations
      depths = obs$depths,
      times = obs$times,
      # label
      name = tag.name
    )
    
    list(dives.processed)
    
  }, depth_files, dive_labels)
  

  names(endpoints) = sapply(dive_labels, function(labels) { labels$tag })
    
  endpoints

}

