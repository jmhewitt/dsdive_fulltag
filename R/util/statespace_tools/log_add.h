//
// Created by Joshua Hewitt on 3/15/21.
//

// Implemented as header-only library for simpler compilation within nimble
// models.

#ifndef ROOKHEADING_H_LOG_ADD_H
#define ROOKHEADING_H_LOG_ADD_H

#include "Rcpp.h"

/**
 * Implement log(1 + e^x) via equation 10 in
 * https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
 */
double log1pexp(double x) {
    if(x <= -37) {
        return exp(x);
    }
    if(x <= 18) {
        return std::log1p(exp(x));
    }
    if(x <= 33.3) {
        return x + exp(-x);
    }
    return x;
}

/*
 * implement log(c) = log(a + b) given log(a) = v, log(b) = data[i],
 * and log(c) = data[i] (updated).  uses the identity:
 *   log(c) = log(a + b) = log(b) + log( 1 + exp(log(a) - log(b)) )
 *
 * corrections are made to account for the possibility that b = 0
 */
double log_add(double log_a, double log_b) {
    if(!std::isfinite(log_a)) {
        if(log_a < 0) {
            return log_b;
        }
    }
    if(!std::isfinite(log_b)) {
        if(log_b < 0) {
            return log_a;
        }
    }
    return log_b + log1pexp(log_a - log_b);
}

// [[Rcpp::export]]
double log_sum_c(std::vector<double> x) {
    auto iter = x.begin();
    auto end = x.end();
    double  res = *(iter++);
    for(iter; iter != end; ++iter)
        res = log_add(res, *iter);
    return res;
}

#endif //ROOKHEADING_H_LOG_ADD_H
