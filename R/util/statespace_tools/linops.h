#ifndef STATESPACE_TOOLS_LINOPS_H
#define STATESPACE_TOOLS_LINOPS_H

// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

#include "log_add.h"

struct NaturalScale {

    /**
     * Componentwise multiplication of \p x by \p scale when \p x and \p scale
     * are stored directly, and overwriting the contents of \p x
     */
    static void cwiseScale(Eigen::VectorXd &x, const Eigen::VectorXd &scale) {
        x.array() *= scale.array();
    }

    /**
     * Post-multiply row-vector \p x by the matrix \p m, i.e., \p x %*% \p m,
     * when \p x and \p m are stored directly, and overwriting the contents of
     * \p x
     * @tparam RowVec Intended to be Eigen::RowVectorXd or
     *   Eigen::Map<Eigen::RowVectorXd>
     */
    template<typename RowVec>
    static void diffuse(RowVec &x, const Eigen::MatrixXd &m) {
        x *= m;
    }

    /**
     * Standardize the vector \p x to have unit 1-Norm when \p x is stored
     * directly, and overwriting the contents of \p x
     */
    static void normalize(Eigen::VectorXd &x) {
        x /= x.sum();
    }

    /**
     * Compute the log-mass of \p x after componentwise reweighting by \p w
     */
    static double logMass(const Eigen::VectorXd &x, const Eigen::VectorXd &w) {
        return std::log((x.array() * w.array()).sum());
    }
};

struct LogScale {

    /**
     * Componentwise multiplication of \p x by \p scale when the log-entries of
     * \p x and \p scale are stored, and overwriting the contents of \p x
     */
    static void cwiseScale(Eigen::VectorXd &x, const Eigen::VectorXd &scale) {
        x += scale;
    }

    static double logMass(const Eigen::VectorXd &x) {
        const double * xiter = x.data();
        const double * xend = xiter + x.size();
        double res = *(xiter++);
        while(xiter != xend) {
            res = log_add(res, *(xiter++));
        }
        return res;
    }

    /**
     * Compute the log-mass of \p x after componentwise reweighting by \p w
     */
    static double logMass(
        const double * xiter,
        const double * witer,
        std::size_t n
    ) {
        const double * xend = xiter + n;
        double res = *(xiter++) + *(witer++);
        while(xiter != xend) {
            res = log_add(res, *(xiter++) + *(witer++));
        }
        return res;
    }

    /**
     * Compute the log-mass of \p x after componentwise reweighting by \p w
     */
    static double logMass(const Eigen::VectorXd &x, const Eigen::VectorXd &w) {
        return logMass(x.data(), w.data(), x.size());
    }

    /**
     * Standardize the vector \p x to have unit 1-Norm when log-entries of \p x
     * are stored, and overwriting the contents of \p x
     */
    static void normalize(Eigen::VectorXd &x) {
        x.array() -= logMass(x);
    }

    /**
     * Post-multiply row-vector \p x by the matrix \p m, i.e., \p x %*% \p m,
     * when the log-entries of \p x and \p m are stored, and overwriting the
     * contents of \p x
     * @tparam RowVec Intended to be Eigen::RowVectorXd or
     *   Eigen::Map<Eigen::RowVectorXd>
     */
    template<typename RowVec>
    static void diffuse(RowVec &x, const Eigen::MatrixXd &m) {

        double * xiter = x.data();
        const double * mcol = m.data();
        std::size_t nrow = x.size();

        RowVec out = x;
        double * oiter = out.data();
        double * oend = oiter + nrow;

        while(oiter != oend) {
            *(oiter++) = logMass(xiter, mcol, nrow);
            mcol += nrow;
        }

        x = out;
    }
};


#endif