/*
* ===========================================================
* File Type: HPP
* File Name: Array2Cube.cpp
* Package Name: RMSS
*
* Created by Anthony-A. Christidis.
* Copyright (c) Anthony-A. Christidis. All rights reserved.
* ===========================================================
*/

// Header files included
#include "Array2Cube.hpp"

arma::cube Array2Cube(Rcpp::NumericVector& my_array) {
  
  Rcpp::NumericVector vec_array(my_array);
  Rcpp::IntegerVector array_dim = vec_array.attr("dim");
  arma::cube cube_array(vec_array.begin(), array_dim[0], array_dim[1], array_dim[2], false);
  return cube_array;
}

arma::ucube Array2UCube(Rcpp::NumericVector& my_array) {
  
  Rcpp::NumericVector vec_array(my_array);
  Rcpp::IntegerVector array_dim = vec_array.attr("dim");
  arma::cube cube_array(vec_array.begin(), array_dim[0], array_dim[1], array_dim[2], false);
  arma::ucube ucube_array = arma::conv_to<arma::ucube>::from(cube_array);
  return ucube_array;
}
