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
  
  Rcpp::IntegerVector array_dim = my_array.attr("dim");
  arma::cube cube_array(my_array.begin(), array_dim[0], array_dim[1], array_dim[2], false);
  return cube_array;
}

arma::ucube Array2UCube(Rcpp::NumericVector& my_array) {
  
  Rcpp::IntegerVector array_dim = my_array.attr("dim");
  arma::ucube ucube_array(array_dim[0], array_dim[1], array_dim[2]);
  double* src_ptr = my_array.begin();
  arma::uword* dst_ptr = ucube_array.memptr();
  arma::uword total_elements = array_dim[0] * array_dim[1] * array_dim[2];
  
  for (arma::uword i = 0; i < total_elements; ++i) {
    dst_ptr[i] = static_cast<arma::uword>(src_ptr[i]);
  }
  
  return ucube_array;
}