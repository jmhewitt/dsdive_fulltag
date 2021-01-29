flatten_tags = function(tag_list, transition_matrices, n_bins, movement_types,
                        pi_discretization, lambda_discretization) {
  
  # initialize flattened structures
  nim_pkg = list(
    data = list(
      depths = NULL,
      times = NULL,
      stages = NULL,
      stage_supports = NULL
    ),
    consts = list(
      segments = NULL,
      transition_matrices = transition_matrices,
      n_bins = n_bins,
      movement_types = movement_types,
      pi_discretization = pi_discretization,
      lambda_discretization = lambda_discretization,
      n_pi = as.integer(pi_discretization[, 'nvals']),
      n_lambda = as.integer(lambda_discretization[, 'nvals']),
      subject_id_labels = NULL
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
    tag_segments = rle(tag$data_gaps)
    valid_segments = which(tag_segments$values == FALSE)
    
    # first indices of observations to analyze
    segment_starts = cumsum(c(1, tag_segments$lengths))
    
    # flatten data
    for(seg_ind in valid_segments) {
      # next available index in nimble package
      flat_ind = length(nim_pkg$data$depths) + 1
      # indices of data segment
      start_ind = segment_starts[seg_ind]
      end_ind = segment_starts[seg_ind + 1] - 1
      seg_inds = seq(from = start_ind, to = end_ind, by = 1)
      # copy data to package
      nim_pkg$data$depths = c(nim_pkg$data$depths, tag$depth.bin[seg_inds])
      nim_pkg$data$stages = c(nim_pkg$data$stages, tag$stages[seg_inds])
      nim_pkg$data$stage_supports = cbind(
        nim_pkg$data$stage_supports,
        tag$stage_support[,seg_inds]
      )
      nim_pkg$data$times = c(nim_pkg$data$times, tag$times[seg_inds])
      # save segment information
      nim_pkg$consts$segments = rbind(
        nim_pkg$consts$segments,
        c(start_ind = flat_ind, length = end_ind - start_ind + 1, 
          subject_id = tag_ind)
      )
    }
    
  }
  
  # compute size constants
  nim_pkg$consts$n_segments = nrow(nim_pkg$consts$segments)
  nim_pkg$consts$n_subjects = length(unique(
    nim_pkg$consts$segments[, 'subject_id']
  ))
  
  nim_pkg
}
