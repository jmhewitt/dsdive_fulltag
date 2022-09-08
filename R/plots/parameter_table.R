parameter_table_script = tar_target(
  name = parameter_table,
  command = {

    model_rep = paste(
      'fit_marginalized_model_', multiple_start_reps, sep =''
    )
    
    # identify posterior parameter samples
    f = dir(
      path = file.path(
        'output', 'mcmc', 'fixed_init_beta', model_rep
      ), 
      full.names = TRUE, 
      pattern = 'mvSamples_[0-9]+'
    )
    
    # get column names for posterior parameter samples
    cnames = readRDS(
      dir(
        path = file.path(
          'output', 'mcmc', 'fixed_init_beta', model_rep
        ), 
        full.names = TRUE, 
        pattern = 'mvSamples_colnames'
      )
    )
    
    # identify parameters to extract
    tgt = c(
      grep(
        pattern = 'lambda', 
        x = cnames, 
        value = TRUE
      ),
      grep(
        pattern = 'pi\\[[13]\\]', 
        x = cnames, 
        value = TRUE
      )
    )
    tgt_inds = match(x = tgt, table = cnames)
    
    # extract samples for target parameters
    samples = mcmc(do.call(rbind, lapply(f, function(f) {
      readRDS(f)[,tgt_inds]
    })))
    colnames(samples) = tgt
    
    # burn = 1:5e3
    burn = 1:7e3

    # check convergence
    for(i in 1:ncol(samples)) {
      print(plot(samples[-burn,i], type = 'l', main = cnames[tgt_inds][i]))
      print(plot(density(samples[-burn,i]), main = cnames[tgt_inds][i]))
    }
    
    
    cbind(
      mean = round(colMeans(mcmc(samples[-burn,])),2),
      sd = round(apply(samples[-burn,],2,sd),3),
      round(HPDinterval(mcmc(samples[-burn,])),2)
    )
    
    effectiveSize(mcmc(samples[-burn,]))
      
    # raw summary table components
    summaries = data.frame(
      param = tgt,
      mean = colMeans(samples),
      sd = apply(samples, 2, sd),
      HPDinterval(samples)
    )
    
    #
    # clean up table
    #
    
    # latex-format parameter names
    summaries$param = paste(
      '$\\',
      gsub(
        pattern = '\\[', 
        replacement = '_', 
        x = gsub(
          pattern = '\\]',
          replacement = '', 
          x = summaries$param
        )
      ),
      '$',
      sep = ''
    )
    
    # format numbers for display
    summaries$mean = format(round(summaries$mean,2))
    summaries$sd = format(round(summaries$sd,2))
    
    # format HPD intervals
    summaries$hpd = paste('(', round(summaries$lower,2), ', ', round(summaries$upper,2), ')', sep = '')
    summaries$lower = NULL
    summaries$upper = NULL
    
    # update column names
    colnames(summaries) = c('Parameter', 'Post. mean', 'Post. s.d.', '95\\% HPD')
    
    # begin printing
    x = print.xtable(
      xtable(
        x = summaries,
        align = 'ccccc',
        caption = 'Posterior distribution summaries for $\\bTheta$ components.', 
        label = 'table:param_ests'
      ),
      booktabs = TRUE,
      sanitize.text.function = identity,
      include.rownames = FALSE,
      caption.placement = 'top',
      comment = FALSE
    )
    
    x = stringr::str_replace_all(
      string = x,
      pattern = '(\\\\bottomrule|\\\\toprule)', 
      replacement = ''
    )
    
    clipr::write_clip(x, allow_non_interactive = TRUE)
    
    # speed differences
    speed_diff_samples = mcmc(apply(samples, 1, function(r) r[2]-r[1]))
    round(mean(speed_diff_samples),2)
    round(HPDinterval(speed_diff_samples),2)
    
    
      0
  }, 
  pattern = map(multiple_start_reps), 
  deployment = 'worker', 
  memory = 'transient',
  storage = 'worker'
)