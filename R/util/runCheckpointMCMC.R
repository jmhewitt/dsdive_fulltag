runCheckpointMCMC = function(
  mcmc, nsamples, out.path, out.name, out.size = Inf, out.time = Inf, 
  verbose = TRUE
) {
  # Parameters:
  #  mcmc - nimble MCMC sampler
  #  nsamples - number of MCMC samples to draw
  #  out.path - directory in which to save sampler output
  #  out.name - base name for output files
  #  out.size - maximum size of sampler output (MB memory) before checkpoint
  #  out.time - ideal maximum time (sec) between checkpoints
  #  verbose - TRUE to output status messages
  
  # initialize number of samples to draw at each checkpoint
  batchSize = nsamples
  
  # validate output directory
  if(!dir.exists(out.path)) {
    stop('Output directory does not exist')
  }
  
  # run sampler once to estimate timings and sizes
  tick = proc.time()[3]
  mcmc$run(niter = 1)
  tock = proc.time()[3]
  
  # extract timing
  sec_sample = tock - tick
  
  # extract initial sample
  samples = as.matrix(mcmc$mvSamples)
  samples2 = as.matrix(mcmc$mvSamples2)
  
  # deconstruct sampler output to estimate sizes
  samples_colnames = colnames(samples)
  samples2_colnames = colnames(samples2)
  colnames(samples) = NULL
  colnames(samples2) = NULL
  
  # extract sizes
  samples_size = lobstr::obj_size(samples)
  samples2_size = lobstr::obj_size(samples2)
  samples_colnames_size = lobstr::obj_size(samples_colnames)
  samples2_colnames_size = lobstr::obj_size(samples2_colnames)
  
  if(verbose) {
    message(paste('Initial sampling rate (sec/sample):', sec_sample))
    message(paste(
      'Estimated sampler duration (sec):', 
      round(sec_sample * nsamples, 2)
    ))
    message(paste(
      'Estimated output size (MB):', 
      round((
        (samples_size + samples2_size) * nsamples + 
        samples_colnames_size + 
        samples2_colnames_size
      ) / 1024^2, 2)
    ))
  }
  
  # batch size constraints due to time and memory constraints
  batchTime = as.numeric(floor(out.time / sec_sample))
  batchMem = as.numeric(floor(
    (out.size * 1024^2 - samples_colnames_size - samples2_colnames_size) / 
    (samples_size + samples2_size)
  ))
  
  # final batch size for checkpoints
  batchSize = min(batchSize, batchTime, batchMem)
  
  # save column headers
  saveRDS(
    samples_colnames, 
    file = file.path(out.path, 
                     paste(out.name, '_mvSamples_colnames.rds', sep = ''))
  )
  saveRDS(
    samples2_colnames, 
    file = file.path(out.path, 
                     paste(out.name, '_mvSamples2_colnames.rds', sep = ''))
  )
  
  # run sampler
  if(verbose) {
    tick = proc.time()[3]
    message(paste(Sys.time(), ' :: Starting Sampling', sep = ''))
  }
  batchInd = 1
  remainingSamples = nsamples
  while(remainingSamples > 0) {
    # determine number of iterations to run
    batch_chunk = min(remainingSamples, batchSize)
    # run sampler
    mcmc$run(niter = batch_chunk, resetMV = TRUE, reset = FALSE)
    # save samples, reducing file size by removing labels
    samples = as.matrix(mcmc$mvSamples)
    samples2 = as.matrix(mcmc$mvSamples2)
    attr(samples, 'dimnames') = NULL
    attr(samples2, 'dimnames') = NULL
    saveRDS(
      samples, 
      file = file.path(out.path, 
                       paste(out.name, '_mvSamples_', batchInd, '.rds', sep = ''))
    )
    saveRDS(
      samples2, 
      file = file.path(out.path, 
                       paste(out.name, '_mvSamples2_', batchInd, '.rds', sep = ''))
    )
    # increment counters
    remainingSamples = remainingSamples - batch_chunk
    batchInd = batchInd + 1
    # output progress
    if(verbose) {
      completed_samples = nsamples - remainingSamples
      message(paste(
        Sys.time(), ' :: Sampling ', 
        round(completed_samples / nsamples * 100, 2),
        '% complete', sep = ''
      ))
      tock = proc.time()[3]
      sec_sample = (tock - tick) / completed_samples
      message(paste(
        '  Estimated time remaining (sec):', 
        round(sec_sample * remainingSamples)
      ))
    }
  }
  
  out.path
}
