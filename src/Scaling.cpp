/*
 * ===========================================================
 * File Type: HPP
 * File Name: Scaling.cpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "Scaling.hpp"

// Functions for the median
double Median(arma::vec& y) {
  return arma::median(y);
}
arma::vec Median(arma::mat& x) {
  return arma::conv_to<arma::vec>::from(arma::median(x));
}
arma::mat MedianData(arma::vec& median_vec, arma::uword& n) {
  
  arma::mat median_x = arma::mat(n, median_vec.size());
  median_x.each_row() = arma::conv_to<arma::rowvec>::from(median_vec);
  return median_x;
}
arma::mat MedianEnsemble(arma::vec& median_vec, arma::uword& n_models) {
  
  arma::mat median_x_ensemble = arma::mat(median_vec.size(), n_models);
  median_x_ensemble.each_col() = median_vec;
  return median_x_ensemble;
}

// Functions for the median absolute deviation (MAD)
double MedianAbsoluteDeviation(arma::vec& y) {
  
  arma::vec y_centered = arma::abs(y - arma::median(y));
  return arma::median(y_centered);
}
arma::vec MedianAbsoluteDeviation(arma::mat& x) {
  
  arma::rowvec median_vec = arma::conv_to<arma::rowvec>::from(arma::median(x));
  arma::mat x_centered = x;
  x_centered.each_row() -= median_vec;
  x_centered = arma::abs(x_centered);
  return arma::conv_to<arma::vec>::from(arma::median(x_centered));;
}
arma::mat MedianAbsoluteDeviationData(arma::vec& mad_vec, arma::uword& n) {
  
  arma::mat mad_x = arma::mat(n, mad_vec.size());
  mad_x.each_row() = arma::conv_to<arma::rowvec>::from(mad_vec);
  return mad_x;
}
arma::mat MedianAbsoluteDeviationEnsemble(arma::vec& mad_vec, arma::uword& n_models) {
  
  arma::mat mad_x_ensemble = arma::mat(mad_vec.size(), n_models);
  mad_x_ensemble.each_col() = mad_vec;
  return mad_x_ensemble;
}

