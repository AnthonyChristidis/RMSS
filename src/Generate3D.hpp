/*
 * ===========================================================
 * File Type: HPP
 * File Name: Generate3D.hpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

#ifndef Generate3D_hpp
#define Generate3D_hpp

// Header files included
#include "EnsembleModel.hpp"

// Function to return 3D vector of intercepts
std::vector<std::vector<std::vector<arma::vec>>> Generate3D_Intercepts(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                       arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                       arma::uword& n_models);

// Function to return 3D vector of coefficients
std::vector<std::vector<std::vector<arma::mat>>> Generate3D_Coefficients(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                         arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                         arma::uword& p, arma::uword& n_models);

// Function to return 3D vector of active samples
std::vector<std::vector<std::vector<arma::umat>>> Generate3D_Active_Samples(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                            arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                            arma::uword& p, arma::uword& n_models);

// Function to return 3D vector of ensemble losses
std::vector<std::vector<std::vector<double>>> Generate3D_Ensemble_Loss(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                       arma::uvec& h, arma::uvec& t, arma::uvec& u);

// Function to initialize 3D vector of CV errors over the folds
std::vector<std::vector<std::vector<arma::vec>>> Generate3D_Prediction_Residuals(arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                                 arma::uword& n);

// Function to return 3D vector of CV error for each tuning parameter configuration
std::vector<std::vector<std::vector<double>>> Generate3D_CV_Error(std::vector<std::vector<std::vector<arma::vec>>> prediction_residuals,
                                                                  arma::uvec& h, arma::uvec& t, arma::uvec& u, 
                                                                  arma::uword& n, arma::uword& n_trim);

#endif // Generate3D_hpp
