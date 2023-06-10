/*
 * ===========================================================
 * File Type: HPP
 * File Name: Scaling.hpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

#ifndef Scaling_hpp
#define Scaling_hpp

// Libraries included
#include <RcppArmadillo.h>

// Functions for the median
double Median(arma::vec& y);
arma::vec Median(arma::mat& x);
arma::mat MedianData(arma::vec& median_vec, arma::uword& n);
arma::mat MedianEnsemble(arma::vec& median_vec, arma::uword& n_models);

// Functions for the median absolute deviation (MAD)
double MedianAbsoluteDeviation(arma::vec& y);
arma::vec MedianAbsoluteDeviation(arma::mat& x);
arma::mat MedianAbsoluteDeviationData(arma::vec& mad_vec, arma::uword& n);
arma::mat MedianAbsoluteDeviationEnsemble(arma::vec& mad_vec, arma::uword& n_models);

#endif // Scaling_hpp
