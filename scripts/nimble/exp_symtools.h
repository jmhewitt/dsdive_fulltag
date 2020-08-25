extern "C"  {

  /*
  Evaluate x = exp(At) * v when provided eigenvectors and eigenvalues for A.

  Parameters:
    evecs - (input) matrix of eigenvectors in column-major format
    evals - (input) vector of eigenvalues
    M - number of rows/cols for evecs
    v - (input) target vector for multiplication
    d - (input) vector of scaling factors
    dInv - (input) vector of inverse scaling factors
    t - value for multiplication
    x - (output) preallocated storage for result of multiplication
    preMultiply - if true, then return v * exp(At) instead of exp(At) * v
*/
void expmAtv_cpp(double* evecs, double* evals, int M, double* v, double* d,
  double* dInv, double t, double* x, bool preMultiply);

/*
  Decompose tridiagonal infinitesimal generator matrix.

  If delta > 0, then the inputs diag, dsuper, and dsub will be modified.

  Parameters:
    diag - (input/output) diagonal from infinitesimal generator on input, and
      eigenvalues on output
    dsuper - (input/output) super-diagonal from infinitesimal generator
    dsub - (input/output) sub-diagonal from infinitesimal generator
    N - dimension of tridiagonal infinitesimal generator matrix
    expm - (output) initial matrix exponential exp(At)
    evecs - (output) storage for eigenvectors in column-major format
    delta - amount of noise to add to sub/super-diagonals to facilitate
      symmetrization
    t - initial value for which to compute exp(At)
*/
void expm_cpp(double* diag, double* dsuper, double* dsub, int N,
  double* expm, double* evecs, double* d, double* dInv, double delta, double t);

}
