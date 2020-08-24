#
# water-drop segmentation function
#

dive.segmentation = function(y, merge.ratio = .5, depth.bins) {
  
  o = order(y, decreasing = TRUE)
  
  #
  # identify and label local modes
  #
  
  # initialize labels for local modes
  modes = numeric(length(y))
  
  # visit datapoints in descending order
  for(ind in o) {
    
    if(modes[ind] == 0) {
      # generate new label for mode
      e = max(modes) + 1
      modes[ind] = e
      
      # flow mode to the left
      if(ind > 1) {
        flow = ind:1
        for(i in 2:length(flow)) {
          if(y[flow[i]] <= y[flow[i-1]]) {
            if(modes[flow[i]] == 0) {
              # associate location with current mode
              modes[flow[i]] = e
            } else {
              # indicate that two modes meet at this location
              modes[flow[i]] = -1
            }
          } else {
            # we have left the basin associated with the current mode
            break
          }
        }
      }
      
      # flow mode to the right
      if(ind < length(y)) {
        flow = ind:length(y)
        for(i in 2:length(flow)) {
          if(y[flow[i]] <= y[flow[i-1]]) {
            if(modes[flow[i]] == 0) {
              # associate location with current mode
              modes[flow[i]] = e
            } else {
              # indicate that two modes meet at this location
              modes[flow[i]] = -1
            }
          } else {
            # we have left the basin associated with the current mode
            break
          }
        }
      }
    }
  }
  
  
  #
  # expand and relabel conflict regions
  #
  
  # relabel conflict locations
  modes[modes==-1] = 0
  
  # identify conflict locations, where merges may potentially occur
  conflict.loci = which(modes==0)
  
  # expand conflict regions
  conflict.label = -1
  if(length(conflict.loci) > 0) {
    for(ind in conflict.loci) {
      # determine all locations where the data attains the conflicted value
      levelset = y == y[ind]
      # uniquely label contiguous regions in the conflict-value set
      levelset.labels = cumsum(c(0, abs(diff(levelset)))) * levelset
      # assign a conflict label to all points contiguous with this conflict
      modes[levelset.labels == levelset.labels[ind]] = conflict.label
      conflict.label = conflict.label - 1
    }
  }
  
  # identify conflict locations, and remove locations for surface observations
  merge.ids = which(modes < 0)
  merge.ids = setdiff(merge.ids, merge.ids[which(y[merge.ids] == 1)])
  nmerge = length(merge.ids)
  
  # merge dives; stop when modes stabilize
  merging = nmerge > 0
  while(merging) {
    
    for(mid in merge.ids) {
      
      # max depth of left-hand segment
      if(mid > 1) {
        # index in (mid:1) of closest mode to mid location
        closest.mode.ind = min(which(modes[mid:1] > 0))
        # label for lhs mode closest to mid
        lhs.mode = modes[mid:1][closest.mode.ind]
        # maximum depth of lhs mode
        lhs.maxdepth = max(y[modes == lhs.mode])
      } else {
        lhs.maxdepth = NA
      }
      
      # max depth of right-hand segment
      if(mid < length(y)) {
        # index in (mid:length(Y)) of closest mode to mid location
        closest.mode.ind = min(which(modes[mid:length(y)] > 0))
        # label for rhs mode closest to mid
        rhs.mode = modes[mid:length(y)][closest.mode.ind]
        # maximum depth of rhs mode
        rhs.maxdepth = max(y[modes == rhs.mode])
      } else {
        rhs.maxdepth = NA
      }
      
      # flag modes to merge if either of the conflict and max depths are similar
      tomerge = max(depth.bins$center[y[mid]] / 
                      depth.bins$center[c(lhs.maxdepth, rhs.maxdepth)], 
                    na.rm = TRUE) > merge.ratio
      
      # merge modes
      if(tomerge) {
        
        # initialize vector of locations to combine
        inds = which(modes == modes[mid])
        
        # initialize label for merged mode
        replacement = NA
        
        # add locations from lhs, if it exists
        if(!is.na(lhs.mode)) {
          inds = c(inds, which(modes == lhs.mode))
          replacement = lhs.mode
        }
        
        # add locations from rhs, if it exists
        if(!is.na(rhs.mode)) {
          inds = c(inds, which(modes == rhs.mode))
          replacement = rhs.mode
        }
        
        # merge modes
        modes[inds] = replacement
        
      }
    }
    
    # identify new conflict locations, and remove locs. for surface obs.
    merge.ids = which(modes < 0)
    merge.ids = setdiff(merge.ids, merge.ids[which(y[merge.ids] == 1)])
    
    
    # stop looping if there are no more conflicts
    if(length(merge.ids) > 0) {
      merging = FALSE
    } else {
      # or if no merges are being made
      if(length(merge.ids) == nmerge) {
        merging = FALSE
      } else {
        nmerge = length(merge.ids)
      }
    }
    
  }
  
  
  #
  # update labels so they are increasing in time
  #
  
  # initialize sequential labels for local modes
  modes.sequential = numeric(length(modes))
  
  # get unique list of non-trivial local modes
  ordered.labels = unique(modes)
  ordered.labels = ordered.labels[ordered.labels > 0]
  
  # assign time-based mode labels
  for(ind in 1:length(ordered.labels)) {
    lab = ordered.labels[ind]
    modes.sequential[modes==lab] = ind
  }
  
  modes.sequential
}
