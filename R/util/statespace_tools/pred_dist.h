#ifndef STATESPACE_TOOLS_PRED_DIST_H
#define STATESPACE_TOOLS_PRED_DIST_H

// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

#include "linops.h"

/**
 * Prediction distribution for state space models provided in an iterator-like
 * container.
 *
 * @tparam LikContainer container providing an iterator that returns
 *   Eigen::VectorXd objects quantifying the sequence of likelihood
 *   distributions [y_i | x_i], where the VectorXd evaluates [y_i | x_i] for
 *   all values of x_i.  As a likelihood, y_i is treated as fixed.
 * @tparam TransitionMatIt container providing an iterator that returns
 *  Eigen::MatrixXd objects quantifying the sequence of transition
 *  matrices [x_{i+1} | x_i]
 * @tparam LinOps implementation of linear algebra operations for scale of data
 */
template<typename LikContainer, typename TxContainer, typename LinOps>
class LatentPrediction {

    /**
     * Recursively compute [x_{i+1} | y_{1:i}] from [x_{i} | y_{1:i-1}].
     * Leaves LikContainer pointing to [y_{i+1} | x_{i+1}] and TransitionMatIt
     * pointing to [x_{i+2} | x_{i+1}].
     */
    void next() {
        // reweight the current prediction distribution by the current obs.
        LinOps::cwiseScale(pred_x, *(lik_it++));
        // forward-diffuse the prediction distribution
        LinOps::diffuse(pred_x_row, *(tx_it++));
        // standardize the distribution
        LinOps::normalize(pred_x);
    }

    // [x_i | y_{1:i-1}]
    Eigen::VectorXd pred_x;

    // row-vector interface to pred_x
    Eigen::Map<Eigen::RowVectorXd> pred_x_row;

    typename LikContainer::iterator lik_it;
    typename TxContainer::iterator tx_it;

public:

    /**
     * @param pred_x0 initial distribution for latent state
     * @param lik
     * @param tx
     */
    LatentPrediction(
        const Eigen::VectorXd & pred_x0, LikContainer & lik,
        TxContainer & tx
    ) : pred_x(pred_x0), pred_x_row(pred_x.data(), pred_x.size()),
        lik_it(lik.begin()), tx_it(tx.begin()) { }

    /**
     * @return reference to current predictive distribution
     */
    Eigen::VectorXd& operator*() { return pred_x; }

    /**
     * prefix increment to advance the predictive distribution
     */
    LatentPrediction& operator++() { next(); return *this; }

};

#endif