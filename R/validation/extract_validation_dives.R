extract_validation_dives = function(tag_list, template_bins, 
                                    validation_pct = 0) {
  
  val_pkg = list(
    consts = list(
      subject_id_labels = NULL
    ),
    data = list(
      dives = NULL
    )
  )
  
  # flatten tag data
  for(tag_ind in 1:length(tag_list)) {
    
    # unwrap tag
    tag = tag_list[[tag_ind]]
    
    # store tag name
    val_pkg$consts$subject_id_labels = c(
      val_pkg$consts$subject_id_labels, tag$tag
    )
    
    # identify all of the observations to analyze
    tag_segments = rle(tag$data_gaps)
    valid_segments = which(tag_segments$values == FALSE)
    
    # first indices of observations to analyze
    segment_starts = cumsum(c(1, tag_segments$lengths))
    
    # last index of baseline period
    baseline_end_ind = max(which(tag$times < tag$baseline_end))
    
    # last pre-exposure index
    last_pre_exposure_ind = max(which(tag$times < tag$exposure_time))
      
    # flatten data
    for(seg_ind in valid_segments) {
      # next available index in nimble package
      flat_ind = length(val_pkg$data$depths) + 1
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
      # if using as a validation dataset, only train using the first portion
      if(validation_pct > 0) {
        train_subset = 1:ceiling(validation_pct * length(seg_inds))
        test_inds = seg_inds[-train_subset]
      }
      # skip segment if there are no transitions to analyze
      if(length(test_inds) == 1) {
        next
      }
      # locally extract data
      depths = tag$depth.bin[test_inds]
      stages = tag$stages[test_inds]
      stage_supports = tag$stage_support[,test_inds]
      times = tag$times[test_inds]
      # build and add covariates
      covariates = rbind(
        intercept = rep(1, length(test_inds)),
        depth =  tag$depths[test_inds],
        deep_depth = tag$depths[test_inds] >= 800,
        shallow_depth = tag$depths[test_inds] < 800,
        non_surface_bin = tag$depth.bin[test_inds] > 1,
        surface_bin = tag$depth.bin[test_inds] == 1,
        time_since_surface = rep(0, length(test_inds)),
        all_shallow_depths_since_surface = rep(0, length(test_inds)),
        daytime = tag$daytime[test_inds],
        moonlit = tag$moonlit[test_inds]
      )
      covariates['time_since_surface',] =  (
        1 + time_in_state(covariates['non_surface_bin',], ncol(covariates))
      ) * covariates['non_surface_bin',]
      all_shallow = TRUE
      for(i in 1:ncol(covariates)) {
        if(covariates['surface_bin',i] == TRUE) {
          all_shallow = TRUE
        } else if(covariates['deep_depth',i] == TRUE) {
          all_shallow = FALSE
        }
        covariates['all_shallow_depths_since_surface',i] = all_shallow
      }
      surface_bin = tag$depth.bin[test_inds] == 1
      
      # use surface observations to identify dive breakpoints
      surface_inds = which(tag$depth.bin[test_inds] == 1)
      # skip segment if there are no clear dives to analyze
      if(length(surface_inds) < 1) {
        next
      }
      # package validation dives
      val_dives = lapply(2:length(surface_inds), function(i) {
        # indices of validation dive
        dive_inds = seq(
          from = surface_inds[i - 1], 
          to = surface_inds[i]
        )
        # nothing interesting to package
        if(length(dive_inds) == 2) {
          return(NULL)
        }
        # package dive
        list(
          depths = depths[dive_inds],
          stages = stages[dive_inds],
          stage_supports = stage_supports[,dive_inds],
          times = times[dive_inds],
          covariates = covariates[,dive_inds],
          surface_bin = surface_bin[dive_inds],
          true_deep = any(covariates['deep_depth', dive_inds] == 1),
          tag = tag$tag,
          dive_start_time = times[dive_inds][1]
        )
      })
      # remove NULL output from validation package
      val_dives = val_dives[!sapply(val_dives, is.null)]
     
      # save segment information
      val_pkg$data$dives = c(val_pkg$data$dives, val_dives)
    }
  }
  
  val_pkg
}
