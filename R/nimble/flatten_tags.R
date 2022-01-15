flatten_tags = function(template_bins, tag_list, depth_threshold, 
                        sattag_timestep, repeated_surface_bin_break = NULL,
                        min_segment_length = 0) {
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
    
    # treat repeated surface bin observations as gaps in animal movement record
    if(!is.null(repeated_surface_bin_break)) {
      # summarize runs of surface vs. non-surface observations
      surface_obs_runs = rle(tag$depth.bin == 1)
      # ids of surface bin runs
      surface_block_ids = which(surface_obs_runs$values == TRUE)
      # identify surface bin runs that are longer than desired
      long_surface_run = surface_obs_runs$lengths[surface_block_ids] >= 
        repeated_surface_bin_break
      # tag interior of long surface runs as data gaps
      if(sum(long_surface_run) > 0) {
        # data indices at which runs begin
        run_starts = cumsum(c(1, surface_obs_runs$lengths))
        # loop over the long surface runs
        for(block_ind in surface_block_ids[long_surface_run]) {
          # label interior of long surface run as a data gap
          tag$gap_after[
            seq(
              from = run_starts[block_ind] + 1,
              by = 1, 
              length.out = surface_obs_runs$lengths[block_ind] - 2
            )
          ] = TRUE
          
        }
      }
    }
    
    # identify all of the observations to analyze
    tag_segments = rle(tag$gap_after)
    valid_segments = which(tag_segments$values == FALSE)
    
    # do not analyze short segments
    valid_segments = valid_segments[
      tag_segments$lengths[valid_segments] > min_segment_length
    ]
    
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
