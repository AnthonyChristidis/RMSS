/*
* ===========================================================
* File Type: HPP
* File Name: Array2Cube.hpp
* Package Name: RMSS
*
* Created by Anthony-A. Christidis.
* Copyright (c) Anthony-A. Christidis. All rights reserved.
* ===========================================================
*/

#ifndef Array2Cube_hpp
#define Array2Cube_hpp

// Libraries included
#include <RcppArmadillo.h>

arma::cube Array2Cube(Rcpp::NumericVector& my_array);
arma::ucube Array2UCube(Rcpp::NumericVector& my_array);

#endif // Array2Cube_hpp
