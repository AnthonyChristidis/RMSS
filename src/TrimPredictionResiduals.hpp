/*
* ===========================================================
* File Type: HPP
* File Name: TrimPredictionResiduals.hpp
* Package Name: RMSS
*
* Created by Anthony-A. Christidis.
* Copyright (c) Anthony-A. Christidis. All rights reserved.
* ===========================================================
*/

#ifndef TrimPredictionResiduals_hpp
#define TrimPredictionResiduals_hpp

// Libraries included
#include <RcppArmadillo.h>

// Function to trim prediction residuals
void TrimPredictionResiduals(arma::vec& prediction_residuals, arma::uvec& sort_order_residuals, arma::uword& n_trim);

#endif // TrimPredictionResiduals_hpp