/*
 * ===========================================================
 * File Type: HPP
 * File Name: SetDifference.hpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

#ifndef SetDifference_hpp
#define SetDifference_hpp

// Libaries included
#include <RcppArmadillo.h>

// Function to return the set difference
arma::uvec SetDifference(const arma::uvec& big, const arma::uvec& small);

#endif // SetDifference_hpp