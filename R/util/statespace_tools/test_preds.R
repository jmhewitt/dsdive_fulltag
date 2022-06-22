if(FALSE) {
  # basic validation that the methods work with transition matrices, etc. 
  # provided directly, and when logged
  
  library(testthat)
  
  M = matrix(c(0,1,1,0),nrow=2)
  # M = matrix(c(0.9,.1,.1,.9),nrow=2)
  
  L = c(.7,.3)
  
  txmats = list(M, M, M, M, M, M, M, M, M, M, M, M, M, M, M, M, M, M, M, M, M)
  liks = list(L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L)
  
  x0 = c(.2,.8)
  
  f = file.path('R', 'util', 'statespace_tools', 'statespace_tools.cpp')
  Rcpp::sourceCpp(f)
  
  expect_equal(
    testPreds(liks = liks, txmats = txmats, x0 = x0),
    lapply(
      testPredsLogEntries(
        liks = lapply(liks, log), 
        txmats = lapply(txmats, log), 
        x0 = log(x0)
      ),
      exp
    )
  )
  
  expect_equal(
    testLL(liks = liks, txmats = txmats, x0 = x0, logEntries = FALSE),
    testLL(
      liks = lapply(liks, log), 
      txmats = lapply(txmats, log), 
      x0 = log(x0), 
      logEntries = TRUE
    )
  )
  
}
