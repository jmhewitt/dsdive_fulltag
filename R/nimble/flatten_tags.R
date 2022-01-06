flatten_tags = function(template_bins, tag_list, depth_threshold, 
                        sattag_timestep) {
  # Parameters:
  #   depth_threshold - depth used to generate prop_recent_deep covariate
  
  # initialize flattened structures
  nim_pkg = list(
    data = list(
      depth_bins = NULL,
      covariates = NULL,
      times = NULL
    ),
    consts = list(
      n_bins = nrow(template_bins),
      widths = 2 * template_bins$halfwidth,
      segments = NULL,
      subject_id_labels = NULL,
      tstep = sattag_timestep
    )
  )
  
  # flatten tag data
  for(tag_ind in 1:length(tag_list)) {
    
    # unwrap tag
    tag = tag_list[[tag_ind]]
    
    # store tag name
    nim_pkg$consts$subject_id_labels = c(
      nim_pkg$consts$subject_id_labels, tag$tag
    )
    
    # identify all of the observations to analyze
    tag_segments = rle(tag$gap_after)
    valid_segments = which(tag_segments$values == FALSE)
    
    # first indices of observations to analyze
    segment_starts = cumsum(c(1, tag_segments$lengths))
    
    # last index of baseline period
    baseline_end_ind = max(which(tag$times < tag$baseline_end))
      
    # flatten data
    for(seg_ind in valid_segments) {
      # next available index in nimble package
      flat_ind = length(nim_pkg$data$depth_bins) + 1
      # start index of data segment
      start_ind = segment_starts[seg_ind]
      # skip segment if it begins after the baseline period
      if(tag$times[start_ind] >= tag$baseline_end) {
        next
      }
      # only analyze baseline portions of tag
      end_ind = min(segment_starts[seg_ind + 1] - 1, baseline_end_ind)
      # indices to analyze
      seg_inds = seq(from = start_ind, to = end_ind, by = 1)
      # skip segment if there are no transitions to analyze
      if(length(seg_inds) <= 1) {
        next
      }
      # copy data to package
      nim_pkg$data$depth_bins = c(nim_pkg$data$depth_bins, 
                                  tag$depth.bin[seg_inds])
      nim_pkg$data$times = c(nim_pkg$data$times, tag$times[seg_inds])
      # add untransformed covariates
      nim_pkg$data$covariates = cbind(
        nim_pkg$data$covariates,
        rbind(
          daytime = tag$daytime[seg_inds],
          moonlit = tag$moonlit[seg_inds],
          time = tag$times[seg_inds],
          depth = tag$depths[seg_inds]
        )
      )
      # save segment information
      nim_pkg$consts$segments = rbind(
        nim_pkg$consts$segments,
        c(start_ind = flat_ind, 
          length = end_ind - start_ind + 1, 
          subject_id = tag_ind,
          end_ind = flat_ind + end_ind - start_ind)
      )
    }
  }
  
  # compute size, etc. constants
  nim_pkg$consts$n_covariates = nrow(nim_pkg$data$covariates)
  nim_pkg$consts$n_segments = nrow(nim_pkg$consts$segments)
  nim_pkg$consts$n_subjects = length(unique(
    nim_pkg$consts$segments[, 'subject_id']
  ))
  
  nim_pkg
}
