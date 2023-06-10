/*
* ===========================================================
* File Type: HPP
* File Name: EnsembleModel.hpp
* Package Name: RMSS
*
* Created by Anthony-A. Christidis.
* Copyright (c) Anthony-A. Christidis. All rights reserved.
* ===========================================================
*/

#ifndef EnsembleModel_hpp
#define EnsembleModel_hpp

// Libraries included
#include <RcppArmadillo.h>
#include <vector>

class EnsembleModel {
  
private:
  
  // Variables supplied by the user
  arma::mat x;
  arma::vec y;
  arma::mat med_x, mad_x;
  arma::mat med_x_ensemble, mad_x_ensemble;
  double med_y, mad_y;
  arma::uword n_models;
  arma::uword h, t, u;
  double tolerance;
  arma::uword max_iter;
  
  // Variables created inside class
  arma::uword n, p;
  arma::mat x_sc; arma::vec y_sc;
  arma::mat coef_mat, coef_mat_candidate;
  double step_size_coef, step_size_trim;
  arma::umat subset_indices, subset_indices_candidate;
  arma::umat active_samples, active_samples_candidate;
  arma::uvec subset_active;
  arma::uvec subset_active_samples;
  arma::uvec group_vec;
  arma::mat final_coef, final_coef_candidate;
  arma::vec final_intercept, final_intercept_candidate;
  arma::vec models_loss, models_loss_candidate;
  double ensemble_loss, ensemble_loss_candidate;
  
  // (+) Functions that update the current state of the ensemble  
  void Compute_Coef(arma::uword& group), Compute_Coef_Candidate(arma::uword& group);
  void Project_Coef(arma::vec& coef_vector, arma::uvec& sort_order_coef);
  void Project_Trim(arma::vec& trim_vector, arma::uvec& sort_order_trim);
  double Compute_Group_Loss(arma::mat& x, arma::vec& y, arma::vec& betas, arma::vec& trim);
  void Update_Subset_Indices(arma::uword& group), Update_Subset_Indices_Candidate(arma::uword& group);
  void Update_Active_Samples(arma::uword& group, arma::vec& new_trim), Update_Active_Samples_Candidate(arma::uword& group, arma::vec& new_trim);
  
public:
  
  // (+) Model Constructor
  EnsembleModel(arma::mat& x, arma::vec& y,
                arma::mat& med_x, arma::mat& mad_x,
                arma::mat& med_x_ensemble, arma::mat& mad_x_ensemble,
                double& med_y, double& mad_y,
                arma::uword& n_models,
                arma::uword& h, arma::uword& t, arma::uword& u,
                double& tolerance,
                arma::uword& max_iter);
  
  // (+) Functions that update the parameters of the ensemble
  void Set_H(arma::uword& h);
  void Set_U(arma::uword& u);
  void Set_T(arma::uword& t);
  void Set_Tolerance(double& tolerance);
  void Set_Max_Iter(arma::uword& max_iter);
  
  // (+) Functions that update the current state of the ensemble  
  void Set_Initial_Indices(arma::umat& subset_indices), Set_Indices_Candidate(arma::umat& subset_indices_candidate);
  void Candidate_Search();
  void Compute_Coef_Ensemble(), Compute_Coef_Ensemble_Candidate();
  void Update_Final_Coef(), Update_Final_Coef_Candidate();
  void Update_Models_Loss(), Update_Models_Loss_Candidate();
  void Update_Ensemble_Loss(), Update_Ensemble_Loss_Candidate();
  void Update_Ensemble();
  
  // (+) Functions that return the current state of the ensemble  
  arma::uvec Get_Model_Subspace(arma::uword& group), Get_Model_Subspace_Candidate(arma::uword& group);
  arma::umat Get_Model_Subspace_Ensemble(), Get_Model_Subspace_Ensemble_Candidate();
  arma::umat Get_Active_Samples();
  arma::vec Get_Final_Intercepts();
  arma::mat Get_Final_Coef();
  double Get_Ensemble_Loss();
  arma::vec Prediction_Residuals_Ensemble(arma::mat& x_test, arma::vec& y_test);
  arma::vec Prediction_Residuals_Models(arma::mat& x_test, arma::vec& y_test);
};

#endif // EnsembleModel_hpp
