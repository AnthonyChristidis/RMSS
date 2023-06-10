/*
* ===========================================================
* File Type: HPP
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

    // Matrix to store candidate indices
    arma::umat candidate_mat = arma::umat(p, n_models);

    // Ensemble losses over all configurations of tuning parameters
    double total_loss_old = 0, total_loss_new = 0;
    for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++)
        for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++)
            for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++)
                total_loss_new += ensembles[h_ind][t_ind][u_ind].Get_Ensemble_Loss();

    // Cycles for neighborhood searches
    do {

        // Setting old total loss
        total_loss_old = total_loss_new;

        //___________________________________________
        // Neighborhood search for non-extreme cases
        //___________________________________________

        // Neighbors for trimming parameter
        for (arma::uword h_ind = 1; h_ind < h.size() - 1; h_ind++) {
            for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
                for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {

                    // LHS Neighbor
                    candidate_mat = ensembles[h_ind - 1][t_ind][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][u_ind].Candidate_Search();

                    // RHS Neighbor
                    candidate_mat = ensembles[h_ind + 1][t_ind][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][u_ind].Candidate_Search();
                }
            }
        }

        // Neighbors of sparsity parameter
        for (arma::uword t_ind = 1; t_ind < t.size() - 1; t_ind++) {
            for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
                for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {

                    // LHS Neighbor
                    candidate_mat = ensembles[h_ind][t_ind - 1][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][u_ind].Candidate_Search();

                    // RHS Neighbor
                    candidate_mat = ensembles[h_ind][t_ind + 1][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][u_ind].Candidate_Search();
                }
            }
        }

        // Neighbors of diversity parameter
        for (arma::uword u_ind = 1; u_ind < u.size() - 1; u_ind++) {
            for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
                for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {

                    // LHS Neighbor
                    candidate_mat = ensembles[h_ind][t_ind][u_ind - 1].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][u_ind].Candidate_Search();

                    // RHS Neighbor
                    candidate_mat = ensembles[h_ind][t_ind][u_ind + 1].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][u_ind].Candidate_Search();
                }
            }
        }

        //_______________________________________
        // Neighborhood search for extreme cases
        //_______________________________________

        // Neighbors for trimming parameter
        if (h.size() > 1) {

            for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {
                for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {

                    // LHS Extreme
                    candidate_mat = ensembles[1][t_ind][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[0][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[0][t_ind][u_ind].Candidate_Search();

                    // RHS Extreme
                    candidate_mat = ensembles[h.size() - 2][t_ind][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[h.size() - 1][t_ind][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h.size() - 1][t_ind][u_ind].Candidate_Search();
                }
            }
        }

        // Neighbors of sparsity parameter
        if (t.size() > 1) {

            for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
                for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++) {

                    // LHS Extreme
                    candidate_mat = ensembles[h_ind][1][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][0][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][0][u_ind].Candidate_Search();

                    // RHS Extreme
                    candidate_mat = ensembles[h_ind][t.size() - 2][u_ind].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t.size() - 1][u_ind].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t.size() - 1][u_ind].Candidate_Search();
                }
            }
        }

        // Neighbors of diversity parameter
        if (u.size() > 1) {

            for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++) {
                for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++) {

                    // LHS Extreme
                    candidate_mat = ensembles[h_ind][t_ind][1].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][0].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][0].Candidate_Search();

                    // RHS Extreme
                    candidate_mat = ensembles[h_ind][t_ind][u.size() - 2].Get_Model_Subspace_Ensemble();
                    ensembles[h_ind][t_ind][u.size() - 1].Set_Indices_Candidate(candidate_mat);
                    ensembles[h_ind][t_ind][u.size() - 1].Candidate_Search();
                }
            }
        }

        // Ensemble losses over all configurations of tuning parameters
        total_loss_new = 0;
        for (arma::uword h_ind = 0; h_ind < h.size(); h_ind++)
            for (arma::uword t_ind = 0; t_ind < t.size(); t_ind++)
                for (arma::uword u_ind = 0; u_ind < u.size(); u_ind++)
                    total_loss_new += ensembles[h_ind][t_ind][u_ind].Get_Ensemble_Loss();

    } while (!(std::abs(total_loss_new - total_loss_old) < neighborhood_search_tolerance));
}

