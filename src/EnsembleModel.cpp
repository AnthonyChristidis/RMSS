/*
 * ===========================================================
 * File Type: CPP
 * File Name: EnsembleModel.cpp
 * Package Name: RMSS
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "EnsembleModel.hpp"

// (+) Model Constructor
EnsembleModel::EnsembleModel(arma::mat& x, arma::vec& y,
                             arma::mat& med_x, arma::mat& mad_x,
                             arma::mat& med_x_ensemble, arma::mat& mad_x_ensemble,
                             double& med_y, double& mad_y,
                             arma::uword& n_models,
                             arma::uword& h, arma::uword& t, arma::uword& u,
                             double& tolerance,
                             arma::uword& max_iter) : 
  x(x), y(y), 
  med_x(med_x), mad_x(mad_x),
  med_x_ensemble(med_x_ensemble), mad_x_ensemble(mad_x_ensemble),
  med_y(med_y), mad_y(mad_y),
  n_models(n_models),
  h(h), t(t), u(u),
  tolerance(tolerance),
  max_iter(max_iter){
  
  // Initialization of dimension of data
  n = x.n_rows;
  p = x.n_cols;
  
  // Initialization of scaled data
  x_sc = (x - med_x) / mad_x;
  y_sc = (y - med_y) / mad_y;
  
  // Initialization of coefficients and indices matrices
  coef_mat = coef_mat_candidate = arma::mat(p, n_models);
  final_intercept = final_intercept_candidate = arma::vec(n_models);
  final_coef = final_coef_candidate = arma::mat(p, n_models);
  subset_indices = subset_indices_candidate = arma::umat(p, n_models);
  active_samples = active_samples_candidate = arma::umat(n, n_models);
  
  // Initialization of vectors for group losses
  models_loss = models_loss_candidate = arma::vec(n_models);
  
  // Initialization of active subsets vectors
  subset_active = arma::uvec(p);
  subset_active_samples = arma::uvec(n);
  
  // Initialization of gradient descent step size for trimming parameter (fixed)
  step_size_trim = 1;
  
  // Initialization of group vector
  group_vec = arma::uvec(1);
}

// (+) Functions that update the parameters of the ensemble
void EnsembleModel::Set_H(arma::uword& h) {
  this->h = h;
}
void EnsembleModel::Set_U(arma::uword& u) {
  this->u = u;
}
void EnsembleModel::Set_T(arma::uword& t) {
  this->t = t;
}
void EnsembleModel::Set_Tolerance(double& tolerance) {
  this->tolerance = tolerance;
}
void EnsembleModel::Set_Max_Iter(arma::uword& max_iter) {
  this->max_iter = max_iter;
}

// (+) Functions that update the current state of the ensemble  
void EnsembleModel::Set_Initial_Indices(arma::umat& subset_indices) {
  this->subset_indices = subset_indices;
}
void EnsembleModel::Set_Indices_Candidate(arma::umat& subset_indices_candidate) {
  this->subset_indices_candidate = subset_indices_candidate;
}
void EnsembleModel::Candidate_Search() {
  
  Compute_Coef_Ensemble_Candidate();
  Update_Final_Coef_Candidate();
  Update_Models_Loss_Candidate();
  Update_Ensemble_Loss_Candidate();
  Update_Ensemble();
}
void EnsembleModel::Compute_Coef_Ensemble() {
  
  for (arma::uword group = 0; group < n_models; group++)
    Compute_Coef(group);
}
void EnsembleModel::Compute_Coef_Ensemble_Candidate() {
  
  for (arma::uword group = 0; group < n_models; group++)
    Compute_Coef_Candidate(group);
}
void EnsembleModel::Compute_Coef(arma::uword& group) {
  
  arma::mat x_subset = x_sc.cols(Get_Model_Subspace(group));
  arma::vec betas, new_betas = arma::zeros(x_subset.n_cols);
  arma::vec trim, new_trim = arma::zeros(n);
  arma::uvec sort_order_coef = arma::uvec(x_subset.n_cols);
  arma::uvec sort_order_trim = arma::uvec(n);
  step_size_coef = 1 /  arma::max(arma::eig_sym(x_subset.t() * x_subset));
  arma::uword iter_count = 0;
  
  do{
    
    // Coefficients update
    betas = new_betas;
    new_betas = betas - step_size_coef * x_subset.t() * (x_subset * betas + new_trim - y_sc);
    Project_Coef(new_betas, sort_order_coef);
    
    // Trimming update
    trim = new_trim;
    new_trim = trim - step_size_trim * (x_subset * betas + new_trim - y_sc);
    Project_Trim(new_trim, sort_order_trim);
    
    // Check if model parmeters have converged
    if (std::abs(Compute_Group_Loss(x_subset, y_sc, new_betas, new_trim) - Compute_Group_Loss(x_subset, y_sc, betas, trim)) < tolerance)
      break;
    
  } while (++iter_count < max_iter);
  
  // Computation of final coefficient based on active sets of predictors and samples
  arma::mat x_final = x_subset.submat(arma::find(new_trim == 0), arma::find(new_betas != 0));
  arma::mat x_final_t_x_final = x_final.t() * x_final;
  new_betas(arma::find(new_betas != 0)) = arma::solve(x_final_t_x_final, arma::eye(arma::size(x_final_t_x_final)), arma::solve_opts::fast) * x_final.t() * y_sc(arma::find(new_trim == 0));
  coef_mat.col(group).zeros();
  group_vec(0) = group;
  coef_mat.submat(Get_Model_Subspace(group), group_vec) = new_betas;
  
  // Updating subset indices for the group
  Update_Subset_Indices(group);
  
  // Updating the active samples for the group
  Update_Active_Samples(group, new_trim);
}
void EnsembleModel::Compute_Coef_Candidate(arma::uword& group) {
  
  arma::mat x_subset = x_sc.cols(Get_Model_Subspace_Candidate(group));
  arma::vec betas, new_betas = arma::zeros(x_subset.n_cols);
  arma::vec trim, new_trim = arma::zeros(n);
  arma::uvec sort_order_coef = arma::uvec(x_subset.n_cols);
  arma::uvec sort_order_trim = arma::uvec(n);
  step_size_coef = 1 / arma::max(arma::eig_sym(x_subset.t() * x_subset));
  arma::uword iter_count = 0;
  
  do {
    
    // Coefficients update
    betas = new_betas;
    new_betas = betas - step_size_coef * x_subset.t() * (x_subset * betas + new_trim - y_sc);
    Project_Coef(new_betas, sort_order_coef);
    
    // Trimming update
    trim = new_trim;
    new_trim = trim - step_size_trim * (x_subset * betas + new_trim - y_sc);
    Project_Trim(new_trim, sort_order_trim);
    
    // Check if model parmeters have converged
    if (std::abs(Compute_Group_Loss(x_subset, y_sc, new_betas, new_trim) - Compute_Group_Loss(x_subset, y_sc, betas, trim)) < tolerance)
      break;
    
  } while (++iter_count < max_iter);
  
  // Computation of final coefficient based on active sets of predictors and samples
  arma::mat x_final = x_subset.submat(arma::find(new_trim == 0), arma::find(new_betas != 0));
  arma::mat x_final_t_x_final = x_final.t() * x_final;
  new_betas(arma::find(new_betas != 0)) = arma::solve(x_final_t_x_final, arma::eye(arma::size(x_final_t_x_final)), arma::solve_opts::fast) * x_final.t() * y_sc(arma::find(new_trim == 0));
  coef_mat_candidate.col(group).zeros();
  group_vec(0) = group;
  coef_mat_candidate.submat(Get_Model_Subspace_Candidate(group), group_vec) = new_betas;
  
  // Updating subset indices for the group
  Update_Subset_Indices_Candidate(group);
  
  // Updating the active samples for the group
  Update_Active_Samples_Candidate(group, new_trim);
}
void EnsembleModel::Project_Coef(arma::vec& coef_vector, arma::uvec& sort_order_coef) {
  
  sort_order_coef = arma::sort_index(arma::abs(coef_vector), "descend");
  coef_vector(sort_order_coef(arma::linspace<arma::uvec>(t, coef_vector.n_elem - 1, coef_vector.n_elem - t))).zeros();
}
void EnsembleModel::Project_Trim(arma::vec& trim_vector, arma::uvec& sort_order_trim) {
  
  sort_order_trim = arma::sort_index(arma::abs(trim_vector), "descend");
  trim_vector(sort_order_trim(arma::linspace<arma::uvec>(n - h, n - 1, h))).zeros();
}
double EnsembleModel::Compute_Group_Loss(arma::mat& x, arma::vec& y, arma::vec& betas, arma::vec& trim) {
  
  return arma::sum(arma::square(y - x * betas - trim));
}
void EnsembleModel::Update_Subset_Indices(arma::uword& group) {
  
  subset_active.zeros();
  subset_active(arma::find(coef_mat.col(group) != 0)).ones();
  subset_indices.col(group) = subset_active;
}
void EnsembleModel::Update_Subset_Indices_Candidate(arma::uword& group) {
  
  subset_active.zeros();
  subset_active(arma::find(coef_mat_candidate.col(group) != 0)).ones();
  subset_indices_candidate.col(group) = subset_active;
}
void EnsembleModel::Update_Active_Samples(arma::uword& group, arma::vec& new_trim) {
  
  subset_active_samples.zeros();
  subset_active_samples(arma::find(new_trim == 0)).ones();
  active_samples.col(group) = subset_active_samples;
}
void EnsembleModel::Update_Active_Samples_Candidate(arma::uword& group, arma::vec& new_trim) {
  
  subset_active_samples.zeros();
  subset_active_samples(arma::find(new_trim == 0)).ones();
  active_samples_candidate.col(group) = subset_active_samples;
}
void EnsembleModel::Update_Final_Coef() {
  
  final_coef = (mad_y * coef_mat) / mad_x_ensemble;
  for(arma::uword group = 0; group < n_models; group++)
    final_intercept(group) =  med_y - arma::as_scalar((final_coef.col(group).t() * med_x_ensemble.col(group)));
}
void EnsembleModel::Update_Final_Coef_Candidate() {
  
  final_coef_candidate = (mad_y * coef_mat_candidate) / mad_x_ensemble;
  for (arma::uword group = 0; group < n_models; group++)
    final_intercept_candidate(group) = med_y - arma::as_scalar((final_coef_candidate.col(group).t() * med_x_ensemble.col(group)));
}
void EnsembleModel::Update_Models_Loss() {
  
  for(arma::uword group = 0; group < n_models; group++)
    models_loss(group) = arma::mean(arma::square(y - final_intercept(group) - x * final_coef.col(group)));
}
void EnsembleModel::Update_Models_Loss_Candidate() {
  
  for (arma::uword group = 0; group < n_models; group++)
    models_loss_candidate(group) = arma::mean(arma::square(y - final_intercept_candidate(group) - x * final_coef_candidate.col(group)));
}
void EnsembleModel::Update_Ensemble_Loss() {
  
  if (u == n_models) {
    
    arma::uword optimal_model = models_loss.index_min();
    for (arma::uword group = 0; group < n_models; group++) {
      
      subset_indices.col(group) = subset_indices.col(optimal_model);
      active_samples.col(group) = active_samples.col(optimal_model);
      final_intercept(group) = final_intercept(optimal_model);
      final_coef.col(group) = final_coef.col(optimal_model);
      models_loss(group) = models_loss(optimal_model);
    }
  }
  
  ensemble_loss = arma::accu(models_loss);
}
void EnsembleModel::Update_Ensemble_Loss_Candidate() {
  
  if (u == n_models) {
    
    arma::uword optimal_model = models_loss.index_min();
    for (arma::uword group = 0; group < n_models; group++) {
      
      subset_indices_candidate.col(group) = subset_indices_candidate.col(optimal_model);
      active_samples_candidate.col(group) = active_samples_candidate.col(optimal_model);
      final_intercept_candidate(group) = final_intercept_candidate(optimal_model);
      final_coef_candidate.col(group) = final_coef_candidate.col(optimal_model);
      models_loss_candidate(group) = models_loss_candidate(optimal_model);
    }
  }
  
  ensemble_loss_candidate = arma::accu(models_loss_candidate);
}
void EnsembleModel::Update_Ensemble() {
  
  if (ensemble_loss_candidate < ensemble_loss) {
    
    subset_indices = subset_indices_candidate;
    active_samples = active_samples_candidate;
    coef_mat = coef_mat_candidate;
    final_intercept = final_intercept_candidate;
    final_coef = final_coef_candidate;
    ensemble_loss = ensemble_loss_candidate;
  }
}

// (+) Functions that return the current state of the ensemble  
arma::uvec EnsembleModel::Get_Model_Subspace(arma::uword& group) {
  return arma::find((arma::sum(subset_indices, 1) - subset_indices.col(group)) < u);
}
arma::uvec EnsembleModel::Get_Model_Subspace_Candidate(arma::uword& group) {
  return arma::find((arma::sum(subset_indices_candidate, 1) - subset_indices_candidate.col(group)) < u);
}
arma::umat EnsembleModel::Get_Model_Subspace_Ensemble() {
  return subset_indices;
}
arma::umat EnsembleModel::Get_Model_Subspace_Ensemble_Candidate() {
  return subset_indices_candidate;
}
arma::umat EnsembleModel::Get_Active_Samples() {
  return active_samples;
}
arma::vec EnsembleModel::Get_Final_Intercepts() {
  return final_intercept;
}
arma::mat EnsembleModel::Get_Final_Coef() {
  return final_coef;
}
double EnsembleModel::Get_Ensemble_Loss() {
  return ensemble_loss;
}
arma::vec EnsembleModel::Prediction_Residuals_Ensemble(arma::mat& x_test, arma::vec& y_test) {
  
  return arma::square(y_test - arma::mean(final_intercept) - arma::mean(x_test * final_coef, 1));
}
arma::vec EnsembleModel::Prediction_Residuals_Models(arma::mat& x_test, arma::vec& y_test) {
  
  arma::vec prediction_residuals = arma::zeros(y_test.n_elem);
  for (arma::uword group = 0; group < n_models; group++) {
    
    prediction_residuals += arma::square(y_test - final_intercept(group) - x_test * final_coef.col(group));
  }
  return prediction_residuals/n_models;
}
