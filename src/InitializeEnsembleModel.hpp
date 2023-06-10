/*
 * ===========================================================
 * File Type: HPP
 * File Name: InitializeEnsembleModel.hpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

#ifndef InitializeEnsembleModel_hpp
#define InitializeEnsembleModel_hpp

// Header files included
#include "EnsembleModel.hpp"

void InitializeEnsembleModel(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                             arma::mat& x, arma::vec& y,
                             arma::mat& med_x, arma::mat& mad_x,
                             arma::mat& med_x_ensemble, arma::mat& mad_x_ensemble,
                             double& med_y, double& mad_y,
                             arma::uword& n_models,
                             arma::uvec& h, arma::uvec& t, arma::uvec& u,
                             double& tolerance,
                             arma::uword& max_iter,
                             arma::umat& initial_split);

#endif // InitializeEnsembleModel_hpp