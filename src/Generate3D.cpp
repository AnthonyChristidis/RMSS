/*
 * ===========================================================
 * File Type: HPP
 * File Name: Generate3D.cpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "Generate3D.hpp"
#include "TrimPredictionResiduals.hpp"

// Function to return 3D vector of intercepts
std::vector<std::vector<std::vector<arma::vec>>> Generate3D_Intercepts(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                       arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                       arma::uword& n_models) {
  
  // Matrix for coefficients
  arma::vec intercepts_vector = arma::vec(n_models);
  
  // Creation of 3D vector for coefficients
  std::vector<std::vector<std::vector<arma::vec>>> intercept3D;
  
  // Allocation of the coefficients
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    
    // 2D vector for a fixed trimming value
    std::vector<std::vector<arma::vec>> intercept3D_h;
    
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      // 1D vector for a fixed sparsity value
      std::vector<arma::vec> intercept3D_t;
      
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
        
        intercepts_vector = ensembles[h_ind][t_ind][u_ind].Get_Final_Intercepts();
        intercept3D_t.push_back(intercepts_vector);
      }
      
      // Adding the 1D vector to the 2D vector of ensembles for a fixed trimming value
      intercept3D_h.push_back(intercept3D_t);
    }
    
    // Adding the 2D vector to the 3D vector of ensembles
    intercept3D.push_back(intercept3D_h);
  }
  
  return intercept3D;
}

// Function to return 3D vector of coefficients
std::vector<std::vector<std::vector<arma::mat>>> Generate3D_Coefficients(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                         arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                         arma::uword& p, arma::uword& n_models) {
  
  // Matrix for coefficients
  arma::mat coef_matrix = arma::mat(p, n_models);
  
  // Creation of 3D vector for coefficients
  std::vector<std::vector<std::vector<arma::mat>>> coef3D;
  
  // Allocation of the coefficients
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    
    // 2D vector for a fixed trimming value
    std::vector<std::vector<arma::mat>> coef3D_h;
    
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      // 1D vector for a fixed sparsity value
      std::vector<arma::mat> coef3D_t;
      
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
        
        coef_matrix = ensembles[h_ind][t_ind][u_ind].Get_Final_Coef();
        coef3D_t.push_back(coef_matrix);
      }
      
      // Adding the 1D vector to the 2D vector of ensembles for a fixed trimming value
      coef3D_h.push_back(coef3D_t);
    }
    
    // Adding the 2D vector to the 3D vector of ensembles
    coef3D.push_back(coef3D_h);
  }
  
  return coef3D;
}

// Function to return 3D vector of active samples
std::vector<std::vector<std::vector<arma::umat>>> Generate3D_Active_Samples(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                            arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                            arma::uword& p, arma::uword& n_models) {
  
  // Matrix for coefficients
  arma::umat coef_matrix = arma::umat(p, n_models);
  
  // Creation of 3D vector for active samples
  std::vector<std::vector<std::vector<arma::umat>>> active3D;
  
  // Allocation of the active samples
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    
    // 2D vector for a fixed trimming value
    std::vector<std::vector<arma::umat>> active3D_h;
    
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      // 1D vector for a fixed sparsity value
      std::vector<arma::umat> active3D_t;
      
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
        
        coef_matrix = ensembles[h_ind][t_ind][u_ind].Get_Active_Samples();
        active3D_t.push_back(coef_matrix);
      }
      
      // Adding the 1D vector to the 2D vector of ensembles for a fixed trimming value
      active3D_h.push_back(active3D_t);
    }
    
    // Adding the 2D vector to the 3D vector of ensembles
    active3D.push_back(active3D_h);
  }
  
  return active3D;
}

// Function to return 3D vector of ensemble losses
std::vector<std::vector<std::vector<double>>> Generate3D_Ensemble_Loss(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                                                                       arma::uvec& h, arma::uvec& t, arma::uvec& u) {
  
  // Creation of 3D vector for ensemble losses
  std::vector<std::vector<std::vector<double>>> loss3D;
  
  // Allocation of the ensemble losses
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    
    // 2D vector for a fixed trimming value
    std::vector<std::vector<double>> loss3D_h;
    
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      // 1D vector for a fixed sparsity value
      std::vector<double> loss3D_t;
      
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
        
        loss3D_t.push_back(ensembles[h_ind][t_ind][u_ind].Get_Ensemble_Loss());
      }
      
      // Adding the 1D vector to the 2D vector of ensembles for a fixed trimming value
      loss3D_h.push_back(loss3D_t);
    }
    
    // Adding the 2D vector to the 3D vector of ensembles
    loss3D.push_back(loss3D_h);
  }
  
  return loss3D;
}

// Function to initialize 3D vector of CV errors over the folds
std::vector<std::vector<std::vector<arma::vec>>> Generate3D_Prediction_Residuals(arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                                 arma::uword& n) {
  
  // Creation of 3D vector for ensemble losses
  std::vector<std::vector<std::vector<arma::vec>>> CVerror3D;
  
  // Allocation of the ensemble losses
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    
    // 2D vector for a fixed trimming value
    std::vector<std::vector<arma::vec>> CVerror3D_h;
    
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      // 1D vector for a fixed sparsity value
      std::vector<arma::vec> CVerror3D_t;
      
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
        
        CVerror3D_t.push_back(arma::zeros(n));
      }
      
      // Adding the 1D vector to the 2D vector of ensembles for a fixed trimming value
      CVerror3D_h.push_back(CVerror3D_t);
    }
    
    // Adding the 2D vector to the 3D vector of ensembles
    CVerror3D.push_back(CVerror3D_h);
  }
  
  return CVerror3D;
}

// Function to return 3D vector of CV error for each tuning parameter configuration
std::vector<std::vector<std::vector<double>>> Generate3D_CV_Error(std::vector<std::vector<std::vector<arma::vec>>> prediction_residuals,
                                                                  arma::uvec& h, arma::uvec& t, arma::uvec& u,
                                                                  arma::uword& n, arma::uword& n_trim) {
  
  // Vectors for trimmed prediction residuals
  arma::vec trimmed_residuals = arma::vec(n);
  arma::uvec sort_order_residuals = arma::uvec(n);
  
  // Creation of 3D vector for ensemble losses
  std::vector<std::vector<std::vector<double>>> CVerror3D;
  
  // Allocation of the ensemble losses
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    
    // 2D vector for a fixed trimming value
    std::vector<std::vector<double>> CVerror3D_h;
    
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      // 1D vector for a fixed sparsity value
      std::vector<double> CVerror3D_t;
      
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
        
        trimmed_residuals = prediction_residuals[h_ind][t_ind][u_ind];
        TrimPredictionResiduals(trimmed_residuals, sort_order_residuals, n_trim);
        CVerror3D_t.push_back(arma::mean(trimmed_residuals));
      }
      
      // Adding the 1D vector to the 2D vector of ensembles for a fixed trimming value
      CVerror3D_h.push_back(CVerror3D_t);
    }
    
    // Adding the 2D vector to the 3D vector of ensembles
    CVerror3D.push_back(CVerror3D_h);
  }
  
  return CVerror3D;
}

