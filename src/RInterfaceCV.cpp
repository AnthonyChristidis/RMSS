/*
 * ===========================================================
 * File Type: CPP
 * File Name: RInterfaceCV.cpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "Scaling.hpp"
#include "Array2Cube.hpp"
#include "SetDifference.hpp"
#include "EnsembleModel.hpp"
#include "InitializeEnsembleModel.hpp"
#include "NeighborhoodSearch.hpp"
#include "Generate3D.hpp"

// [[Rcpp::export]]
Rcpp::List RInterfaceCV(arma::mat& x, arma::vec& y,
                        arma::uword& n_models,
                        arma::uvec& h, arma::uvec& t, arma::uvec& u,
                        double& tolerance,
                        arma::uword& max_iter,
                        Rcpp::NumericVector& initial_split_array,
                        arma::umat& initial_split,
                        arma::uword& neighborhood_search,
                        double& neighborhood_search_tolerance,
                        arma::uword& n_folds,
                        double& alpha,
                        double& gamma,
                        arma::uword& n_threads) {
  
  // Storing number of samples and predictors
  arma::uword n = x.n_rows;
  arma::uword p = x.n_cols;
  
  // Storing number of trimmed samples for CV procedure
  arma::uword n_trim = std::round(n * alpha);
  
  //____________________________
  // Cross-validation procedure
  //____________________________
  
  // Transformation of array of initial splits to cube
  arma::ucube initial_split_folds = Array2UCube(initial_split_array);
  
  // Creating indices for the folds of the data
  const arma::uvec indin = arma::linspace<arma::uvec>(0, n - 1, n);
  const arma::uvec inint = arma::linspace<arma::uvec>(0, n, n_folds + 1);
  
  // Creating 3D vector to store prediction residuals of CV procedure
  std::vector<std::vector<std::vector<arma::vec>>> prediction_residuals = Generate3D_Prediction_Residuals(h, t, u, n);
  
  // Vector for trimming grid of the fold
  arma::uvec h_fold = arma::uvec(h.n_elem);
  
  // Multithreading over the folds of the data
  // # pragma omp parallel for num_threads(n_threads)
  for (arma::uword fold = 0; fold < n_folds; fold++) {
    
    // Creating the training and test data
    arma::uvec test = arma::linspace<arma::uvec>(inint[fold], inint[fold + 1] - 1, inint[fold + 1] - inint[fold]);
    arma::uvec train = SetDifference(indin, test);
    arma::mat x_train = x.rows(train);
    arma::mat x_test = x.rows(test);
    arma::vec y_train = y.elem(train);
    arma::vec y_test = y.elem(test);
    
    // Scaling data
    arma::uword n_in_fold = x_train.n_rows;
    arma::vec med_x = Median(x_train);
    arma::mat med_x_data = MedianData(med_x, n_in_fold);
    arma::mat med_x_ensemble = MedianEnsemble(med_x, n_models);
    double med_y = Median(y_train);
    arma::vec mad_x = MedianAbsoluteDeviation(x_train);
    arma::mat mad_x_data = MedianAbsoluteDeviationData(mad_x, n_in_fold);
    arma::mat mad_x_ensemble = MedianAbsoluteDeviationEnsemble(mad_x, n_models);
    double mad_y = MedianAbsoluteDeviation(y_train);
    
    // Creation of 3D vector to store ensembles
    std::vector<std::vector<std::vector<EnsembleModel>>> ensembles;
    
    // Trimming grid for the fold
    for(arma::uword h_id = 0; h_id < h_fold.n_elem; h_id++)
      h_fold[h_id] = std::round((double)n_in_fold / n * h[h_id]);

    // Initialization of ensembles
    InitializeEnsembleModel(ensembles,
                            x_train, y_train,
                            med_x_data, mad_x_data,
                            med_x_ensemble, mad_x_ensemble,
                            med_y, mad_y,
                            n_models,
                            h_fold, t, u,
                            tolerance,
                            max_iter,
                            initial_split_folds.slice(fold));
    
    // Neighborhood search
    if (neighborhood_search)
      NeighborhoodSearch(ensembles,
                         h, t, u,
                         p, n_models,
                         neighborhood_search_tolerance);
    
    // Storing the prediction residuals of the CV procedure
    for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++)
      for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++)
        for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) 
          prediction_residuals[h_ind][t_ind][u_ind](test) = gamma * ensembles[h_ind][t_ind][u_ind].Prediction_Residuals_Ensemble(x_test, y_test) +
            (1 - gamma) * ensembles[h_ind][t_ind][u_ind].Prediction_Residuals_Models(x_test, y_test);
  }
  
  //_____________________________
  // Ensemble model on full data
  //_____________________________
  
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
  
  //____________________
  // CV Function output
  //____________________
  
  // Creating list for output
  Rcpp::List output;
  output["cv_error"] = Generate3D_CV_Error(prediction_residuals, h, t, u, n, n_trim);
  output["active_samples"] = Generate3D_Active_Samples(ensembles, h, t, u, p, n_models);
  output["intercepts"] = Generate3D_Intercepts(ensembles, h, t, u, n_models);
  output["coef"] = Generate3D_Coefficients(ensembles, h, t, u, p, n_models);
  output["loss"] = Generate3D_Ensemble_Loss(ensembles, h, t, u);
  
  // Return final output
  return output;
}
