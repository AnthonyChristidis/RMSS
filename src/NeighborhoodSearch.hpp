/*
 * ===========================================================
 * File Type: HPP
 * File Name: NeighborhoodSearch.hpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

#ifndef NeighborhoodSearch_hpp
#define NeighborhoodSearch_hpp

// Libraries included
#include <RcppArmadillo.h>
#include <vector>

// Header files included
#include "EnsembleModel.hpp"

void NeighborhoodSearch(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                        arma::uvec& h, arma::uvec& t, arma::uvec& u,
                        arma::uword& p, arma::uword& n_models,
                        double& neighborhood_search_tolerance);

#endif // NeighborhoodSearch_hpp