#
# compile package
#

base_dir = file.path(getwd(), 'R', 'util', 'expokit', 'rexpokit_bare')

system(paste(
  'R CMD COMPILE',
  paste('CXXFLAGS=', '"',
        paste('-I', 
              file.path(find.package(c('Rcpp')), 'include'), 
              sep = '', collapse = ' '), 
        '"', sep = ''),
  'FFLAGS=-fallow-argument-mismatch',
  file.path(base_dir, 'expokit_flattened.f'),
  sep = ' '
))

expokit_gpadm = nimble::nimbleExternalCall(
  prototype = function(ideg = integer(1), m = integer(1), t = double(1), 
                       H = double(2), ldh = integer(1), wsp = double(1), 
                       lwsp = integer(1), ipiv = integer(1), iexph = integer(1),
                       ns = integer(1), iflag = integer(1)){},
  returnType = void(),
  Cfun = 'wrapdgpadm_',
  headerFile = file.path(base_dir, 'expokit.h'), 
  oFile = file.path(base_dir, 'expokit_flattened.o')
)

expocall_gpadm = nimble::nimbleFunction(
  run = function(H = double(2), t = double(0), nrows = integer(0), 
                 ncols = integer(0)) {
    
    returnType(double(2))
    
    ideg <- integer(6, length = 1)
    m <- integer(nrows, length = 1)
    ldh <- integer(ncols, length = 1)
    lwsp <- integer(4*m*m+ideg+1 + 1, length = 1)
    wsp <- numeric(0, length = lwsp[1])
    ipiv <- integer(0, length = m[1])
    iexph <- integer(1)
    ns <- integer(1)
    iflag <- integer(1)
    timpl <- numeric(t, length = 1)
    
    expokit_gpadm(ideg, m, timpl, H, ldh, wsp, lwsp, ipiv, iexph, ns, iflag)
    
    ind_start <- iexph[1]
    ind_end <- ind_start + nrows * nrows - 1

    expm <- matrix(wsp[ind_start:ind_end], nrow = nrows, ncol = ncols)

    return(expm)
  }
)