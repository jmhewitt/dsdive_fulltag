#include "statespace_tools.h"

#include "pred_dist.h"
#include "marginal_lik.h"
#include "dictionary_decoder.h"

// [[Rcpp::export]]
std::vector<Eigen::VectorXd> testPreds(
    std::vector<Eigen::VectorXd> liks,
    std::vector<Eigen::MatrixXd> txmats,
    Eigen::VectorXd x0
) {

    LatentPrediction<
        std::vector<Eigen::VectorXd>,
        std::vector<Eigen::MatrixXd>
    > pred(x0, liks, txmats);

    std::vector<Eigen::VectorXd> pred_vecs;

    for(int i = 0; i < liks.size(); ++i)
        pred_vecs.push_back(*(++pred));

    return pred_vecs;
}

// [[Rcpp::export]]
double testLL(
    std::vector<Eigen::VectorXd> liks,
    std::vector<Eigen::MatrixXd> txmats,
    Eigen::VectorXd x0
) {
    return ll_marginal(x0, liks, txmats, liks.size());
}

/**
 * Marginal likelihood for a single-layer discrete-time, discrete-space
 * state space model.  The number of states is assumed to be small enough that
 * all numerical operations can be implemented using dense matrix methods.
 *
 * Furthermore, it is assumed that there is considerable redundancy in both the
 * state transition and state observation distributions.  The unique
 * distributions are passed in via "dictionary" objects to reduce memory
 * requirements for datasets with long observation records.  The dictionary
 * compression also helps reduce the computational cost for generating the
 * distributions for all timepoints, but this step is done before this method
 * is used.
 *
 * @param obs_lik_dict m \times n matrix where each column contains the
 *   observation distribution [y | x] for a fixed y.  The format implies there
 *   are m latent states for x = 0,...,m-1 and n unique values for
 *   y = 0,...,n-1.
 * @param obs sequence of observed values for y
 * @param txmat_dict data vector for an m \times m \times q array that stores
 *   q unique state transition matrices for the state space model.
 * @param txmat_seq Sequence of transition matrices used in the state space
 *   model.
 * @param x0 initial distribution of latent states
 * @return
 */
// [[Rcpp::export]]
double llCompressed(
    Eigen::MatrixXd obs_lik_dict, std::vector<int> obs,
    std::vector<double> txmat_dict, std::vector<int> txmat_seq,
    Eigen::VectorXd x0
) {

    unsigned int nstates = x0.size();

    // decoder for the transition matrices
    MatrixMapper txmap(txmat_dict.data(), nstates, nstates);
    DictionaryDecoder<
        MatrixMapper,
        Eigen::Map<Eigen::MatrixXd>,
        int
        > txmats(txmap, txmat_seq);

    // decoder for the likelihood functions
    ColumnMapper likmap(obs_lik_dict);
    DictionaryDecoder<
        ColumnMapper,
        Eigen::Map<Eigen::VectorXd>,
        int
    > liks(likmap, obs);

    return ll_marginal(x0, liks, txmats, obs.size());
}

double nimLLCompressedRaw(
    double* obs_lik_dict, double* obs, double* txmat_dict, double* txmat_seq,
    double* x0, int nstates, int nt
) {

    // decoder for the transition matrices
    MatrixMapper txmap(txmat_dict, nstates, nstates);
    DictionaryDecoder<
        MatrixMapper,
        Eigen::Map<Eigen::MatrixXd>,
        double
    > txmats(txmap, txmat_seq, nt);

    // decoder for the likelihood functions
    ColumnMapper likmap(obs_lik_dict, nstates);
    DictionaryDecoder<
        ColumnMapper,
        Eigen::Map<Eigen::VectorXd>,
        double
    > liks(likmap, obs, nt);

    // wrap prior for initial state
    Eigen::Map<Eigen::VectorXd> x0vec(x0, nstates);

    return ll_marginal(x0vec, liks, txmats, nt);
}

/**
 * Marginal likelihood for a two-layer discrete-time, discrete-space
 * state space model.  The number of states is assumed to be small enough that
 * all numerical operations can be implemented using dense matrix methods.
 *
 * Furthermore, it is assumed that there is considerable redundancy in both the
 * state transition and state observation distributions.  The unique
 * distributions are passed in via "dictionary" objects to reduce memory
 * requirements for datasets with long observation records.  The dictionary
 * compression also helps reduce the computational cost for generating the
 * distributions for all timepoints, but this step is done before this method
 * is used.
 *
 * @param obs_lik_dict An array of m n \times n transition matrices where
 *   each n \times n transition matrix contains the transition distribution
 *   [y_{i+1} | y_i, x] for a fixed y_i and y_{i+1}.  The format implies there
 *   are m latent states for x = 0,...,m-1 and n unique values for
 *   y = 0,...,n-1.  The array is assumed to be stored in column major format
 *   with indices specifying (x, y_i, y_{i+1}).  The order of the indices is
 *   chosen s.t. a consecutive sequence of m elements can be used as the
 *   likelihood [y_{i+1} | y_i, \cdot] in filtering computations.
 * @param obs sequence of observed values for y
 * @param txmat_dict data vector for an m \times m \times q array that stores
 *   q unique state transition matrices for the state space model, providing
 *   [x_{i+1} | x_i].
 * @param txmat_seq Sequence of transition matrices used in the state space
 *   model.
 * @param x0 initial distribution of latent states
 * @return
 */
// [[Rcpp::export]]
double ll2LayerCompressed(
    std::vector<double> obs_lik_dict, std::vector<int> obs,
    std::vector<double> txmat_dict, std::vector<int> txmat_seq,
    Eigen::VectorXd x0, int num_obs_states
) {

    unsigned int nstates = x0.size();

    // decoder for the hidden state's transition matrices
    MatrixMapper txmap(txmat_dict.data(), nstates, nstates);
    DictionaryDecoder<
        MatrixMapper,
        Eigen::Map<Eigen::MatrixXd>,
        int
    > txmats(txmap, txmat_seq);

    // decoder for the observed state's transition matrices
    ColumnMapper3 likmap(obs_lik_dict.data(), nstates, num_obs_states);
    DictionaryDecoder<
        ColumnMapper3,
        Eigen::Map<Eigen::VectorXd>,
        int,
        MarkovDictIter<ColumnMapper3, Eigen::Map<Eigen::VectorXd>, int>
    > liks(likmap, obs);

    return ll_marginal(x0, liks, txmats, obs.size() - 1);
}

double nimLL2LayerCompressedRaw(
    double* obs_lik_dict, double* obs,
    double* txmat_dict, double* txmat_seq,
    double* x0, int num_obs_states, int num_latent_states, int nt
) {

    unsigned int nstates = num_latent_states;

    // decoder for the hidden state's transition matrices
    MatrixMapper txmap(txmat_dict, nstates, nstates);
    DictionaryDecoder<
      MatrixMapper,
      Eigen::Map<Eigen::MatrixXd>,
      double
    > txmats(txmap, txmat_seq, nt);

    // decoder for the observed state's transition matrices
    ColumnMapper3 likmap(obs_lik_dict, nstates, num_obs_states);
    DictionaryDecoder<
      ColumnMapper3,
      Eigen::Map<Eigen::VectorXd>,
      double,
    MarkovDictIter<ColumnMapper3, Eigen::Map<Eigen::VectorXd>, double>
    > liks(likmap, obs, nt);

    // wrap prior for initial state
    Eigen::Map<Eigen::VectorXd> x0vec(x0, nstates);
    
    return ll_marginal(x0vec, liks, txmats, nt - 1);
}

/**
 * Final (marginal) prediction distribution [x_n | y_{1:n}] for a two-layer
 * discrete-time, discrete-space state space model.  The number of states is
 * assumed to be small enough that all numerical operations can be implemented
 * using dense matrix methods.
 *
 * Furthermore, it is assumed that there is considerable redundancy in both the
 * state transition and state observation distributions.  The unique
 * distributions are passed in via "dictionary" objects to reduce memory
 * requirements for datasets with long observation records.  The dictionary
 * compression also helps reduce the computational cost for generating the
 * distributions for all timepoints, but this step is done before this method
 * is used.
 *
 * @param obs_lik_dict An array of m n \times n transition matrices where
 *   each n \times n transition matrix contains the transition distribution
 *   [y_{i+1} | y_i, x] for a fixed y_i and y_{i+1}.  The format implies there
 *   are m latent states for x = 0,...,m-1 and n unique values for
 *   y = 0,...,n-1.  The array is assumed to be stored in column major format
 *   with indices specifying (x, y_i, y_{i+1}).  The order of the indices is
 *   chosen s.t. a consecutive sequence of m elements can be used as the
 *   likelihood [y_{i+1} | y_i, \cdot] in filtering computations.
 * @param obs sequence of observed values for y
 * @param txmat_dict data vector for an m \times m \times q array that stores
 *   q unique state transition matrices for the state space model, providing
 *   [x_{i+1} | x_i].
 * @param txmat_seq Sequence of transition matrices used in the state space
 *   model.
 * @param x0 initial distribution of latent states
 * @return
 */
// [[Rcpp::export]]
Eigen::VectorXd finalPred2LayerCompressed(
    std::vector<double> obs_lik_dict, std::vector<int> obs,
    std::vector<double> txmat_dict, std::vector<int> txmat_seq,
    Eigen::VectorXd x0, int num_obs_states
) {

    unsigned int nstates = x0.size();

    typedef DictionaryDecoder<
                MatrixMapper,
                Eigen::Map<Eigen::MatrixXd>,
                int
            > TxContainer;

    typedef DictionaryDecoder<
                ColumnMapper3,
                Eigen::Map<Eigen::VectorXd>,
                int,
                MarkovDictIter<ColumnMapper3, Eigen::Map<Eigen::VectorXd>, int>
            > LikContainer;

    // decoder for the hidden state's transition matrices
    MatrixMapper txmap(txmat_dict.data(), nstates, nstates);
    TxContainer txmats(txmap, txmat_seq);

    // decoder for the observed state's transition matrices
    ColumnMapper3 likmap(obs_lik_dict.data(), nstates, num_obs_states);
    LikContainer liks(likmap, obs);

    // initialize container and updating methods for predictive distribution
    LatentPrediction<LikContainer, TxContainer> pred_dist(x0, liks, txmats);

    // iterate to compute the predictive distribution for the final state
    unsigned int nobs = obs.size() - 1;
    for(unsigned int i = 0; i < nobs; ++i) {
        ++pred_dist;
    }

    // extract [x_n | y_{1:n-1}], which becomes [x_n | y_{1:n}] by assuming
    // a flat likelihood at time y_n, i.e., that there is no information on the
    // final observed transition
    Eigen::VectorXd xf = *pred_dist;

    return xf;
}
