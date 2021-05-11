eda_targets = list(
  
  tar_target(
    name = mle_speeds, 
    command = {
      
      # unwrap tag
      tag = raw_sattags[[1]]
      
      # depth bin widths and count
      n_bins = nrow(template_bins)
      bin_widths = 2 * template_bins$halfwidth
      
      # 
      nobs = 1e2
      
      tick = proc.time()[3]
      mle_params = do.call(rbind, lapply(1:nobs, function(ind) {
        
        # time between observations (sec)
        tstep = diff(as.numeric(tag$times[ind + 0:1]))
        
        # depth bin at start and end of observation
        s0 = tag$depth.bin[ind]
        sf = tag$depth.bin[ind + 1]
        
        # # optimize transition parameters
        # o = optim(c(0,0), fn = function(theta) {
        #   expm(tstep * buildInfinitesimalGenerator(
        #     pi = plogis(theta[1]), 
        #     lambda = exp(theta[2]), 
        #     M = n_bins, 
        #     stage = 6, # unrestricted bin transitions
        #     widths = bin_widths
        #   ))[s0, sf]
        # }, control = list(fnscale = -1, maxit = 1e3), method = 'BFGS')
        
        piseq = c(.01, .5, .99)
        
        # profile-optimization of transition parameters
        o = lapply(piseq, function(pi) {
          optim(0,  fn = function(log_lambda) {
            expm(tstep * buildInfinitesimalGenerator(
              pi = pi, 
              lambda = exp(log_lambda), 
              M = n_bins, 
              stage = 6, # unrestricted bin transitions
              widths = bin_widths
            ))[s0, sf]
          }, control = list(fnscale = -1), method = 'BFGS')
        })
        
        if(any(sapply(o, function(o) o$convergence) != 0)) {
          warning(paste('Convergence failed for obs_ind', ind))
        }
        
        o_ind = which.max(sapply(o, function(o) o$value))
        
        # extract and back-transform parameters
        data.frame(pi = piseq[o_ind], lambda = exp(o[[o_ind]]$par))
      }))
      tock = proc.time()[3]
      
      tock-tick
      
      mle_params$t = tag$times[1:nobs]
      mle_params$ind = 1:nobs
      mle_params$depth = tag$depths[1:nobs]
      mle_params$depth.next = tag$depths[1:nobs + 1]
      mle_params$apparent_speed = abs(mle_params$depth.next - mle_params$depth)/300
      mle_params$depth.bin = tag$depth.bin[1:nobs]
      mle_params$depth.bin.next = tag$depth.bin[1:nobs + 1]
      mle_params$depth.bin.prev = c(NA, tag$depth.bin[1:(nobs-1)])
      mle_params$depth.bin.txd = abs(mle_params$depth.bin.next - mle_params$depth.bin)
        
      ggplot(mle_params %>% pivot_longer(cols = c(pi,lambda, depth), 
                                         names_to = 'param'), 
             aes(x = t, y = value)) + 
        geom_line() + 
        geom_point() + 
        theme_few() + 
        facet_wrap(~param, scales = 'free', ncol = 1)
      
      
      library(viridis)
      
      ggplot(mle_params %>% 
               mutate(pi = factor(pi),
                      lambda = cut(lambda, breaks = c(-.1,0.5,1,6))) %>%
               pivot_longer(cols = c(pi,lambda), names_to = 'parameter'), 
             aes(x = t, y = depth, 
                 group = 1,
                 col = value)) +
                 # col = factor(pi))) +
                 # col = cut(lambda, breaks = 5))) + 
        facet_wrap(~parameter, ncol = 1) + 
        geom_line() + 
        geom_point() + 
        scale_y_reverse() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2') +
        theme_few() 
      
      
      ggplot(mle_params, aes(x = depth.bin.next, y = lambda)) + 
        geom_point() + 
        scale_x_continuous(breaks = 1:16) 
      
      ggplot(mle_params, aes(x = factor(depth.bin), y = lambda,
                             col = depth.bin >= depth.bin.prev)) + 
        geom_boxplot() +
        geom_point() 
      
      ggplot(mle_params, aes(x = factor(pi), y = lambda)) + 
        geom_boxplot() +
        geom_point() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2')
      
      ggplot(mle_params, aes(x = factor(pi), y = lambda),
                             col = depth.bin >= depth.bin.prev) + 
        geom_boxplot() +
        geom_point() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2')
      
      
      browser()
    }, 
    pattern = map(raw_sattags)
  )
  
)
