/*
 * ===========================================================
 * File Type: CPP
 * File Name: RInterface.cpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "Scaling.hpp"
#include "EnsembleModel.hpp"
#include "InitializeEnsembleModel.hpp"
#include "NeighborhoodSearch.hpp"
#include "Generate3D.hpp"

// [[Rcpp::export]]
Rcpp::List RInterface(arma::mat& x, arma::vec& y,
                      arma::uword& n_models,
                      arma::uvec& h, arma::uvec& t, arma::uvec& u,
                      double& tolerance,
                      arma::uword& max_iter,
                      arma::umat& initial_split,
                      arma::uword& neighborhood_search,
                      double& neighborhood_search_tolerance) {
  
  // Storing number of samples and predictors
  arma::uword n = x.n_rows;
  arma::uword p = x.n_cols;
  
  // Scaling data
  arma::vec med_x = Median(x);
  arma::mat med_x_data = MedianData(med_x, n);
  arma::mat med_x_ensemble = MedianEnsemble(med_x, n_models);
  double med_y = Median(y);
  arma::vec mad_x = MedianAbsoluteDeviation(x);
  arma::mat mad_x_data = MedianAbsoluteDeviationData(mad_x, n);
  arma::mat mad_x_ensemble = MedianAbsoluteDeviationEnsemble(mad_x, n_models);
  double mad_y = MedianAbsoluteDeviation(y);
  
  // Creation of 3D vector to store ensembles
  std::vector<std::vector<std::vector<EnsembleModel>>> ensembles;
  
  // Initialization of ensembles
  InitializeEnsembleModel(ensembles,
                          x, y,
                          med_x_data, mad_x_data,
                          med_x_ensemble, mad_x_ensemble,
                          med_y, mad_y,
                          n_models,
                          h, t, u,
                          tolerance,
                          max_iter,
                          initial_split);
  
  // Neighborhood search
  if (neighborhood_search)
    NeighborhoodSearch(ensembles,
                       h, t, u,
                       p, n_models,
                       neighborhood_search_tolerance);
  
  // Creating list for output
  Rcpp::List output;
  output["active_samples"] = Generate3D_Active_Samples(ensembles, h, t, u, p, n_models);
  output["intercepts"] = Generate3D_Intercepts(ensembles, h, t, u, n_models);
  output["coef"] = Generate3D_Coefficients(ensembles, h, t, u, p, n_models);
  output["loss"] = Generate3D_Ensemble_Loss(ensembles, h, t, u);
  
  // Return final output
  return output;
}
