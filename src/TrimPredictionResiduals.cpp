/*
* ===========================================================
* File Type: HPP
* File Name: TrimPredictionResiduals.cpp
* Package Name: RMSS
*
* Created by Anthony-A. Christidis.
* Copyright (c) Anthony-A. Christidis. All rights reserved.
* ===========================================================
*/

// Header files included
#include "TrimPredictionResiduals.hpp"

// Function to trim prediction residuals
void TrimPredictionResiduals(arma::vec& prediction_residuals, arma::uvec& sort_order_residuals, arma::uword& n_trim) {

    sort_order_residuals = arma::sort_index(arma::abs(prediction_residuals), "descend");
    prediction_residuals(sort_order_residuals(arma::linspace<arma::uvec>(0, n_trim - 1, n_trim))).zeros();
}

