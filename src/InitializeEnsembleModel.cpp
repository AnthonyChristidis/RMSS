/*
 * ===========================================================
 * File Type: HPP
 * File Name: InitializeEnsembleModel.cpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "InitializeEnsembleModel.hpp"

void InitializeEnsembleModel(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                             arma::mat& x, arma::vec& y,
                             arma::mat& med_x, arma::mat& mad_x,
                             arma::mat& med_x_ensemble, arma::mat& mad_x_ensemble,
                             double& med_y, double& mad_y,
                             arma::uword& n_models,
                             arma::uvec& h, arma::uvec& t, arma::uvec& u,
                             double& tolerance,
                             arma::uword& max_iter,
                             arma::umat& initial_split) {
  
  // 3D initialization of the ensembles
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    
    // 2D vector for a fixed trimming value
    std::vector<std::vector<EnsembleModel>> ensembles_h;
    
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      // 1D vector for a fixed sparsity value
      std::vector<EnsembleModel> ensembles_t;
      
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++)
        ensembles_t.push_back(EnsembleModel(x, y,
                                            med_x, mad_x,
                                            med_x_ensemble, mad_x_ensemble,
                                            med_y, mad_y,
                                            n_models,
                                            h(h_ind), t(t_ind), u(u_ind),
                                            tolerance,
                                            max_iter));
      
      // Adding the 1D vector to the 2D vector of ensembles for a fixed trimming value
      ensembles_h.push_back(ensembles_t);
    }
    
    // Adding the 2D vector to the 3D vector of ensembles
    ensembles.push_back(ensembles_h);
  }
  
  // Adding the initial split for the case where there is no overlap between the models
  // Also compute the coefficients
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      
      ensembles[h_ind][t_ind][0].Set_Initial_Indices(initial_split);
      ensembles[h_ind][t_ind][0].Compute_Coef_Ensemble();
      ensembles[h_ind][t_ind][0].Update_Final_Coef();
      ensembles[h_ind][t_ind][0].Update_Models_Loss();
      ensembles[h_ind][t_ind][0].Update_Ensemble_Loss();
    }
  }
  
  // Compute ensemble coefficients for cases with overlap
  arma::umat candidate_mat = arma::umat(x.n_cols, n_models);
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      for (arma::uword u_ind = 1; u_ind < u.size(); u_ind++) {
        
        candidate_mat = ensembles[h_ind][t_ind][u_ind - 1].Get_Model_Subspace_Ensemble();
        ensembles[h_ind][t_ind][u_ind].Set_Initial_Indices(candidate_mat);
        ensembles[h_ind][t_ind][u_ind].Compute_Coef_Ensemble();
        ensembles[h_ind][t_ind][u_ind].Update_Final_Coef();
        ensembles[h_ind][t_ind][u_ind].Update_Models_Loss();
        ensembles[h_ind][t_ind][u_ind].Update_Ensemble_Loss();
      }
    }
  }
}

