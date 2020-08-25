#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "exp_symtools.h"

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

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
  double* dInv, double t, double* x, bool preMultiply = false) {

    Map<MatrixXd> evecs_mat(evecs, M, M);
    Map<VectorXd> evals_vec(evals, M);

    MatrixXd unscaled = evecs_mat *
      (evals_vec.array() * t).exp().matrix().asDiagonal() *
      evecs_mat.transpose();

    // TODO: roll in the better correction

    // correct small values
    double* raw = unscaled.data();
    double n = unscaled.rows() * unscaled.cols();
    for(int i=0; i<n; i++) {
      double tmp = *raw;
      *(raw++) = std::abs(tmp);
    }

    Map<VectorXd> res(x, M);
    Map<VectorXd> v_vec(v, M);
    Map<VectorXd> d_vec(d, M);
    Map<VectorXd> dInv_vec(dInv, M);

    if(preMultiply) {
      res = v_vec.transpose() * d_vec.asDiagonal() * unscaled *
        dInv_vec.asDiagonal();
    } else {
      res = d_vec.asDiagonal() * unscaled * dInv_vec.asDiagonal() * v_vec;
    }

}

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
  double* expm, double* evecs, double* d, double* dInv, double delta = 1e-15,
  double t = 300) {

  // VectorXd diag = A.diagonal();
  // VectorXd dsuper = A.diagonal(1);
  // VectorXd dsub = A.diagonal(-1);

  // make A symmetrizable by adding small noise
  if(delta > 0) {

    // shift main diagonal
    double delta2 = 2 * delta;
    for(int i=0; i < N; i++)
      diag[i] -= delta2;

    // shift super and sub diagonals
    for(int i=0; i<N-1; i++) {
      dsuper[i] += delta;
      dsub[i] += delta;
    }
    dsuper[0] += delta;
    dsub[N-2] += delta;

  }

  // direct construction of similarity-transformed symmetric version of matrix
  VectorXd E(N-1);
  for(int i=0; i<N-1; i++) {
    E(i) = sqrt(dsuper[i] * dsub[i]);
  }

  // LAPACK internals
  char compz = 'I';
  Map<MatrixXd> evecs_mat(evecs, N, N);
  VectorXd work(2*N-2);
  int info = 0;

  F77_CALL(dsteqr)(&compz, &N, diag, E.data(), evecs, &N, work.data(), &info);

   // build similarity transform
   Map<VectorXd> d_vec(d, N);
   Map<VectorXd> dInv_vec(dInv, N);
   d[0] = 1.0;
   dInv[0] = 1.0;
   for(int i=1; i<N; i++) {
     double tmp = std::sqrt(dsub[i-1] / dsuper[i-1]) * d[i-1];
     d[i] = tmp;
     dInv[i] = 1/tmp;
   }

   Map<VectorXd> diag_vec(diag, N);

   MatrixXd unscaled = evecs_mat *
     (diag_vec.array() * t).exp().matrix().asDiagonal() *
     evecs_mat.transpose();

   // "correct" small, negative values
   double* raw = unscaled.data();
   int entries = N*N;
   for(int i=0; i<entries; i++) {
     double tmp = *raw;
     *(raw++) = std::abs(tmp);
   }

   Map<MatrixXd> expm_mat(expm, N, N);
   expm_mat = d_vec.asDiagonal() * unscaled * dInv_vec.asDiagonal();

}
