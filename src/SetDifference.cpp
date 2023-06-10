/*
* ===========================================================
* File Type: HPP
* File Name: SetDifference.cpp
* Package Name: RMSS
*
* Created by Anthony-A. Christidis.
* Copyright (c) Anthony-A. Christidis. All rights reserved.
* ===========================================================
*/

// Header files included
#include "SetDifference.hpp"

// Function to return the set difference
arma::uvec SetDifference(const arma::uvec& big, const arma::uvec& small) {

    // Find set difference between a big and a small set of variables.
    // Note: small is a subset of big (both are sorted).
    int m = small.n_elem;
    int n = big.n_elem;
    arma::uvec test = arma::uvec(n, arma::fill::zeros);
    arma::uvec zeros = arma::uvec(n - m, arma::fill::zeros);

    for (int j = 0; j < m; j++) {
        test[small[j]] = small[j];
    }

    test = big - test;
    if (small[0] != 0)
        test[0] = 1;

    zeros = arma::find(test != 0);
    return(zeros);
}

