#ifndef STATESPACE_TOOLS_MARGINAL_LIK_H
#define STATESPACE_TOOLS_MARGINAL_LIK_H

#include "pred_dist.h"

/**
 * Evaluate the marginal likelihood for a statespace model
 *
 * @tparam LikContainer container providing an iterator that returns
 *   Eigen::VectorXd objects quantifying the sequence of likelihood
 *   distributions [y_i | x_i], where the VectorXd evaluates [y_i | x_i] for
 *   all values of x_i.  As a likelihood, y_i is treated as fixed.  Note that
 *   the definition of y_i and x_i depend on the calling function.  For example,
 *   it is possible that y_i may actually represent a pair of observations, if
 *   one if fitting a model that has a Markov evolution for "y_i" as well.
 * @tparam TransitionMatIt container providing an iterator that returns
 *  Eigen::MatrixXd objects quantifying the sequence of transition
 *  matrices [x_{i+1} | x_i]
 */
template<typename LikContainer, typename TxContainer>
double ll_marginal(
    const Eigen::VectorXd & pred_x0, LikContainer & lik, TxContainer & tx,
    unsigned int nobs
) {

    double ll = 0;

    LatentPrediction<LikContainer, TxContainer> pred_dist(pred_x0, lik, tx);

    auto lik_it = lik.begin();
    for(unsigned int i = 0; i < nobs; ++i) {
        ll += std::log(((*pred_dist).array() * (*lik_it).array()).sum());
        // mathematically, must update likelihood after prediction distribution
        ++pred_dist;
        ++lik_it;
    }

    return ll;
}

#endif