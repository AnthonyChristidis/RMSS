/*
 * ===========================================================
 * File Type: CPP
 * File Name: NeighborhoodSearch.cpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "NeighborhoodSearch.hpp"

void NeighborhoodSearch(std::vector<std::vector<std::vector<EnsembleModel>>>& ensembles,
                        arma::uvec& h, arma::uvec& t, arma::uvec& u,
                        arma::uword& p, arma::uword& n_models,
                        double& neighborhood_search_tolerance) {
  
  // Pre-compute initial total loss
  double total_loss_old = 0;
  std::vector<std::vector<std::vector<double>>> cached_losses(h.size(),
                                                              std::vector<std::vector<double>>(t.size(), std::vector<double>(u.size())));
  
  for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
    for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
      for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
        cached_losses[h_ind][t_ind][u_ind] = ensembles[h_ind][t_ind][u_ind].Get_Ensemble_Loss();
        total_loss_old += cached_losses[h_ind][t_ind][u_ind];
      }
    }
  }
  
  // Track which configurations changed this iteration
  std::vector<std::vector<std::vector<bool>>> changed(h.size(),
                                                      std::vector<std::vector<bool>>(t.size(), std::vector<bool>(u.size(), false)));
  
  // Neighborhood search iterations
  do {
    double total_loss_new = 0;
    std::fill(changed.begin(), changed.end(), 
              std::vector<std::vector<bool>>(t.size(), std::vector<bool>(u.size(), false)));
    
    // Helper lambda to try a neighbor and track changes
    auto try_neighbor = [&](arma::uword h_curr, arma::uword t_curr, arma::uword u_curr,
                            arma::uword h_neighbor, arma::uword t_neighbor, arma::uword u_neighbor) {
      double old_loss = ensembles[h_curr][t_curr][u_curr].Get_Ensemble_Loss();
      
      // Try the neighbor configuration
      arma::umat candidate_mat = ensembles[h_neighbor][t_neighbor][u_neighbor].Get_Model_Subspace_Ensemble();
      ensembles[h_curr][t_curr][u_curr].Set_Indices_Candidate(candidate_mat);
      ensembles[h_curr][t_curr][u_curr].Candidate_Search();
      
      double new_loss = ensembles[h_curr][t_curr][u_curr].Get_Ensemble_Loss();
      
      // Track if this configuration improved
      if (std::abs(new_loss - old_loss) > 1e-10) {
        changed[h_curr][t_curr][u_curr] = true;
        cached_losses[h_curr][t_curr][u_curr] = new_loss;
      }
    };
    
    //___________________________________________
    // Consolidated neighborhood search
    //___________________________________________
    
    for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
      for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
        for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
          
          // Check h-dimension neighbors
          if (h_ind > 0) {
            try_neighbor(h_ind, t_ind, u_ind, h_ind - 1, t_ind, u_ind);
          }
          if (h_ind < h.size() - 1) {
            try_neighbor(h_ind, t_ind, u_ind, h_ind + 1, t_ind, u_ind);
          }
          
          // Check t-dimension neighbors  
          if (t_ind > 0) {
            try_neighbor(h_ind, t_ind, u_ind, h_ind, t_ind - 1, u_ind);
          }
          if (t_ind < t.size() - 1) {
            try_neighbor(h_ind, t_ind, u_ind, h_ind, t_ind + 1, u_ind);
          }
          
          // Check u-dimension neighbors
          if (u_ind > 0) {
            try_neighbor(h_ind, t_ind, u_ind, h_ind, t_ind, u_ind - 1);
          }
          if (u_ind < u.size() - 1) {
            try_neighbor(h_ind, t_ind, u_ind, h_ind, t_ind, u_ind + 1);
          }
        }
      }
    }
    
    // Compute new total loss using cached values
    for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
      for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
        for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {
          total_loss_new += cached_losses[h_ind][t_ind][u_ind];
        }
      }
    }
    
    // Check for convergence
    if (std::abs(total_loss_new - total_loss_old) < neighborhood_search_tolerance) {
      break;
    }
    
    total_loss_old = total_loss_new;
    
  } while (true);
}