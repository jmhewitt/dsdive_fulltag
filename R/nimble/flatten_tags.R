flatten_tag_data = function(depth_files, dive_endpoints, template_bins, tag_sex,
                            cee_starts, validation_proportion) {
  
  # initialize flattened structures
  nim_pkg = list(
    data = list(
      depths = NULL,
      times = NULL
    ),
    consts = list(
      endpoint_priors = NULL,
      dive_relations = NULL,
      dive_relations_validation = NULL,
      tag_covariates = NULL
    )
  )
  
  # merge tag records
  for(i in 1:length(dive_endpoints)) {
    
    attach(dive_endpoints[[i]])
    
    # extract sex of tagged whale
    nim_pkg$consts$tag_covariates = c(
      nim_pkg$consts$tag_covariates,
      as.numeric(tag_sex %>% dplyr::filter(deployid == name) %>% 
                   dplyr::select(sex)) - 1
    )
    
    # determine prior distribution's temporal support for tag
    tag_start = endpoint.inds$t_lwr[dive.ranges$T0_endpoint][1] 
    tag_end = endpoint.inds$t_upr[dive.ranges$T3_endpoint]
    tag_end = tag_end[length(tag_end)]
    
    # determine which CEE's overlap with tag record
    cees_experienced = which((tag_start <= cee_starts) & (cee_starts <= tag_end))
    
    
    #
    # id's of dives to keep from record
    #
    
    # keep deep dives and dives that pass earlier data completeness tests
    dive.ids = base::intersect(which(dive.flags),
                         dive.ranges$dive.id[dive.ranges$type == 'Deep'])
    
    # keep dives with more than three observations
    dive.ids = base::intersect(dive.ids, 
                         which(dive.ranges$end.ind - dive.ranges$start.ind + 1 > 
                                 3))
    
    # keep dives that end before an exposure event
    if(length(cees_experienced) > 0) {
      exclude_after_time = cee_starts[min(cees_experienced)]
      dive.ids = base::intersect(dive.ids,
                           which(endpoint.inds$t_lwr[dive.ranges$T3_endpoint] <= 
                                   exclude_after_time))
    }
    
    # keep dives with observations that start and end at the surface
    dive.ids = base::intersect(
      dive.ids,
      which(depths[endpoint.inds$ind[dive.ranges$T0_endpoint]] == 1 & 
            depths[endpoint.inds$ind[dive.ranges$T3_endpoint]] == 1)
    )
    
    # make train/test partitions
    if(validation_proportion > 0) {
      
      # partition dives by id
      test_dives = sample(x = dive.ids, 
                          size = length(dive.ids) * validation_proportion)
      dive.ids = base::setdiff(dive.ids, test_dives)
      
    } else {
      test_dives = NULL
    }
    
    # id's of endpoints to keep from record
    endpoint.ids = endpoint.inds %>% 
      dplyr::filter(dive_end %in% dive.ids | dive_start %in% dive.ids) %>% 
      dplyr::select(endpoint.id) %>% 
      unlist()
    
    # initialize simple storage and lookup for merged id's
    endpoint.inds$merged.id = NA
    
    # merge endpoints
    for(j in 1:nrow(endpoint.inds)) {
      # only process valid endpoints
      if(endpoint.inds$endpoint.id[j] %in% endpoint.ids) {
        
        # construct id for endpoint in merged dataset
        flat_endpoint_id = ifelse(is.null(nim_pkg$consts$endpoint_priors), 0,
                                  nrow(nim_pkg$consts$endpoint_priors)) + 1
        
        # store merged endpoint id
        endpoint.inds$merged.id[j] = flat_endpoint_id
        
        # merge endpoint
        nim_pkg$consts$endpoint_priors = rbind(
          nim_pkg$consts$endpoint_priors,
          c(t_lwr = endpoint.inds$t_lwr[j], 
            t_upr = endpoint.inds$t_upr[j])
        )
      }
    }
    
    # merge each dive
    for(j in 1:nrow(dive.ranges)) {
      # only process valid dives
      if(dive.ranges$dive.id[j] %in% c(dive.ids, test_dives)) {
        
        # get index at which first depth will be stored
        flat_depth_ind = length(nim_pkg$data$depths) + 1
        
        # merge observed depths
        nim_pkg$data$depths = c(
          nim_pkg$data$depths,
          depths[dive.ranges$start.ind[j]:dive.ranges$end.ind[j]]
        )
        
        # merge observation times
        nim_pkg$data$times = c(
          nim_pkg$data$times,
          times[dive.ranges$start.ind[j]:dive.ranges$end.ind[j]]
        )
        
        # merge dive record
        if(dive.ranges$dive.id[j] %in% dive.ids) { 
          
          # construct id for dive in merged dataset
          flat_dive_id = ifelse(is.null(nim_pkg$consts$dive_relations), 0, 
                                nrow(nim_pkg$consts$dive_relations)) + 1
          
          # merge dive record
          nim_pkg$consts$dive_relations = rbind(
            nim_pkg$consts$dive_relations,
            c(T0_endpoint = endpoint.inds$merged.id[dive.ranges$T0_endpoint[j]], 
              T3_endpoint = endpoint.inds$merged.id[dive.ranges$T3_endpoint[j]],
              tag = i,
              depth_first = flat_depth_ind,
              depth_last = length(nim_pkg$data$depths))
          )
          
        } else {
          
          # construct id for dive in merged validation dataset
          flat_dive_id = ifelse(
            is.null(nim_pkg$consts$dive_relations_validation), 0, 
            nrow(nim_pkg$consts$dive_relations_validation)
          ) + 1
          
          # merge dive record
          nim_pkg$consts$dive_relations_validation = rbind(
            nim_pkg$consts$dive_relations_validation,
            c(tag = i,
              depth_first = flat_depth_ind,
              depth_last = length(nim_pkg$data$depths))
          )
          
        }
        
      }
    }
    
    detach(dive_endpoints[[i]])
    rm(endpoint.inds)
    
  }
  
  # locations where there are consecutive deep dives
  length(which(
    nim_pkg$consts$dive_relations[-1,'T0_endpoint'] == 
      nim_pkg$consts$dive_relations[1:(nrow(nim_pkg$consts$dive_relations)-1),
                                    'T3_endpoint']
  ))
  
  # extract sizes
  nim_pkg$consts$N_tags = length(dive_endpoints)
  nim_pkg$consts$N_dives = nrow(nim_pkg$consts$dive_relations)
  nim_pkg$consts$N_endpoints = nrow(nim_pkg$consts$endpoint_priors)
  nim_pkg$consts$N_bins = nrow(template_bins)
  
  # export bin widths
  nim_pkg$consts$widths = template_bins$halfwidth * 2
  
  if(validation_proportion==0) {
    nim_pkg$consts$dive_relations_validation = NULL
  }
  
  nim_pkg
}
