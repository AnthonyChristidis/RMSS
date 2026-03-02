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
#include <algorithm>

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
  
  // NEW: Initialize cache variables
  cached_step_sizes.resize(n_models, -1.0);
  cached_subspaces.resize(n_models);
  subspace_cache_valid.resize(n_models, false);
  row_sums_subset = arma::zeros<arma::uvec>(p);
  row_sums_subset_candidate = arma::zeros<arma::uvec>(p);
  row_sums_cache_valid = false;
  row_sums_candidate_cache_valid = false;
}

// NEW: Cache management functions
void EnsembleModel::invalidate_cache() {
  row_sums_cache_valid = false;
  std::fill(subspace_cache_valid.begin(), subspace_cache_valid.end(), false);
  std::fill(cached_step_sizes.begin(), cached_step_sizes.end(), -1.0);
}

void EnsembleModel::invalidate_candidate_cache() {
  row_sums_candidate_cache_valid = false;
}

void EnsembleModel::update_row_sums_cache() {
  if (!row_sums_cache_valid) {
    row_sums_subset = arma::sum(subset_indices, 1);
    row_sums_cache_valid = true;
  }
}

void EnsembleModel::update_row_sums_candidate_cache() {
  if (!row_sums_candidate_cache_valid) {
    row_sums_subset_candidate = arma::sum(subset_indices_candidate, 1);
    row_sums_candidate_cache_valid = true;
  }
}

// (+) Functions that update the parameters of the ensemble
void EnsembleModel::Set_H(arma::uword& h) {
  this->h = h;
}
void EnsembleModel::Set_U(arma::uword& u) {
  this->u = u;
  invalidate_cache();
  invalidate_candidate_cache();
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
  invalidate_cache();
}
void EnsembleModel::Set_Indices_Candidate(arma::umat& subset_indices_candidate) {
  this->subset_indices_candidate = subset_indices_candidate;
  invalidate_candidate_cache();
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

// OPTIMIZED: Main coefficient computation function
void EnsembleModel::Compute_Coef(arma::uword& group) {
  
  arma::uvec model_subspace = Get_Model_Subspace(group);
  arma::mat x_subset = x_sc.cols(model_subspace);
  
  // Cache step size computation - only compute once per subspace
  if (cached_step_sizes[group] < 0 || !subspace_cache_valid[group]) {
    arma::mat XtX = x_subset.t() * x_subset;
    
    // ROBUST STEP SIZE COMPUTATION WITH MULTIPLE FALLBACKS
    bool step_size_computed = false;
    
    // Method 1: Try original eig_sym approach
    try {
      arma::vec eigenvals = arma::eig_sym(XtX);
      double max_eig = arma::max(eigenvals);
      if (max_eig > 1e-12 && std::isfinite(max_eig)) {
        cached_step_sizes[group] = 1.0 / max_eig;
        step_size_computed = true;
      }
    } catch (...) {
      // eig_sym failed, continue to next method
    }
    
    // Method 2: Try with regularization if original failed
    if (!step_size_computed) {
      try {
        arma::mat XtX_reg = XtX;
        XtX_reg.diag() += 1e-8;
        arma::vec eigenvals = arma::eig_sym(XtX_reg);
        double max_eig = arma::max(eigenvals);
        if (max_eig > 1e-12 && std::isfinite(max_eig)) {
          cached_step_sizes[group] = 1.0 / max_eig;
          step_size_computed = true;
        }
      } catch (...) {
        // Regularized eig_sym also failed, continue to next method
      }
    }
    
    // Method 3: Try SVD approach
    if (!step_size_computed) {
      try {
        arma::vec singular_values;
        bool svd_success = arma::svd(singular_values, XtX);
        if (svd_success && singular_values.n_elem > 0 && 
            singular_values(0) > 1e-12 && std::isfinite(singular_values(0))) {
          cached_step_sizes[group] = 1.0 / singular_values(0);
          step_size_computed = true;
        }
      } catch (...) {
        // SVD also failed, continue to fallback
      }
    }
    
    // Method 4: Trace-based fallback
    if (!step_size_computed) {
      try {
        double trace_val = arma::trace(XtX);
        if (trace_val > 1e-12 && std::isfinite(trace_val) && XtX.n_rows > 0) {
          double avg_eigenval = trace_val / XtX.n_rows;
          cached_step_sizes[group] = 1.0 / avg_eigenval;
          step_size_computed = true;
        }
      } catch (...) {
        // Even trace failed
      }
    }
    
    // Method 5: Ultimate conservative fallback
    if (!step_size_computed) {
      cached_step_sizes[group] = 1e-15;  // Very conservative step size
    }
    
    subspace_cache_valid[group] = true;
  }
  step_size_coef = cached_step_sizes[group];
  
  arma::vec betas = arma::zeros(x_subset.n_cols);
  arma::vec new_betas = betas;
  arma::vec trim = arma::zeros(n);
  arma::vec new_trim = trim;
  
  // Pre-compute frequently used quantities
  arma::vec Xty = x_subset.t() * y_sc;
  arma::mat XtX = x_subset.t() * x_subset;
  
  arma::uword iter_count = 0;
  double prev_loss = std::numeric_limits<double>::max();
  
  do {
    betas = new_betas;
    trim = new_trim;
    
    // More efficient gradient computation
    arma::vec residuals = x_subset * betas + trim - y_sc;
    
    // Coefficients update - reuse pre-computed matrices
    new_betas = betas - step_size_coef * (XtX * betas + x_subset.t() * trim - Xty);
    Project_Coef(new_betas);
    
    // Trimming update
    new_trim = trim - step_size_trim * residuals;
    Project_Trim(new_trim);
    
    // More efficient convergence check
    double current_loss = arma::dot(residuals, residuals);
    if (std::abs(current_loss - prev_loss) < tolerance) break;
    prev_loss = current_loss;
    
  } while (++iter_count < max_iter);
  
  // Optimized final coefficient computation
  arma::uvec active_predictors = arma::find(new_betas != 0);
  arma::uvec active_samples = arma::find(new_trim == 0);
  
  if (!active_predictors.empty() && !active_samples.empty()) {
    arma::mat X_active = x_subset.submat(active_samples, active_predictors);
    arma::vec y_active = y_sc.elem(active_samples);
    
    // Use normal equations with regularization for numerical stability
    arma::mat XtX_active = X_active.t() * X_active;
    arma::vec Xty_active = X_active.t() * y_active;
    XtX_active.diag() += 1e-8;  // Small regularization
    
    arma::vec beta_active = arma::solve(XtX_active, Xty_active, arma::solve_opts::fast);
    new_betas.elem(active_predictors) = beta_active;
  }
  
  coef_mat.col(group).zeros();
  group_vec(0) = group;
  coef_mat.submat(model_subspace, group_vec) = new_betas;
  
  // Updating subset indices for the group
  Update_Subset_Indices(group);
  
  // Updating the active samples for the group
  Update_Active_Samples(group, new_trim);
}

// OPTIMIZED: Candidate coefficient computation
void EnsembleModel::Compute_Coef_Candidate(arma::uword& group) {
  
  arma::uvec model_subspace = Get_Model_Subspace_Candidate(group);
  arma::mat x_subset = x_sc.cols(model_subspace);
  
  // ROBUST STEP SIZE COMPUTATION WITH MULTIPLE FALLBACKS
  arma::mat XtX = x_subset.t() * x_subset;
  bool step_size_computed = false;
  
  // Method 1: Try original eig_sym approach
  try {
    arma::vec eigenvals = arma::eig_sym(XtX);
    double max_eig = arma::max(eigenvals);
    if (max_eig > 1e-12 && std::isfinite(max_eig)) {
      step_size_coef = 1.0 / max_eig;
      step_size_computed = true;
    }
  } catch (...) {
    // eig_sym failed, continue to next method
  }
  
  // Method 2: Try with regularization if original failed
  if (!step_size_computed) {
    try {
      arma::mat XtX_reg = XtX;
      XtX_reg.diag() += 1e-8;
      arma::vec eigenvals = arma::eig_sym(XtX_reg);
      double max_eig = arma::max(eigenvals);
      if (max_eig > 1e-12 && std::isfinite(max_eig)) {
        step_size_coef = 1.0 / max_eig;
        step_size_computed = true;
      }
    } catch (...) {
      // Regularized eig_sym also failed, continue to next method
    }
  }
  
  // Method 3: Try SVD approach
  if (!step_size_computed) {
    try {
      arma::vec singular_values;
      bool svd_success = arma::svd(singular_values, XtX);
      if (svd_success && singular_values.n_elem > 0 && 
          singular_values(0) > 1e-12 && std::isfinite(singular_values(0))) {
        step_size_coef = 1.0 / singular_values(0);
        step_size_computed = true;
      }
    } catch (...) {
      // SVD also failed, continue to fallback
    }
  }
  
  // Method 4: Trace-based fallback
  if (!step_size_computed) {
    try {
      double trace_val = arma::trace(XtX);
      if (trace_val > 1e-12 && std::isfinite(trace_val) && XtX.n_rows > 0) {
        double avg_eigenval = trace_val / XtX.n_rows;
        step_size_coef = 1.0 / avg_eigenval;
        step_size_computed = true;
      }
    } catch (...) {
      // Even trace failed
    }
  }
  
  // Method 5: Ultimate conservative fallback
  if (!step_size_computed) {
    step_size_coef = 0.01;  // Very conservative step size
  }
  
  arma::vec betas = arma::zeros(x_subset.n_cols);
  arma::vec new_betas = betas;
  arma::vec trim = arma::zeros(n);
  arma::vec new_trim = trim;
  
  // Pre-compute frequently used quantities
  arma::vec Xty = x_subset.t() * y_sc;
  
  arma::uword iter_count = 0;
  double prev_loss = std::numeric_limits<double>::max();
  
  do {
    betas = new_betas;
    trim = new_trim;
    
    // More efficient gradient computation
    arma::vec residuals = x_subset * betas + trim - y_sc;
    
    // Coefficients update
    new_betas = betas - step_size_coef * (XtX * betas + x_subset.t() * trim - Xty);
    Project_Coef(new_betas);
    
    // Trimming update
    new_trim = trim - step_size_trim * residuals;
    Project_Trim(new_trim);
    
    // More efficient convergence check
    double current_loss = arma::dot(residuals, residuals);
    if (std::abs(current_loss - prev_loss) < tolerance) break;
    prev_loss = current_loss;
    
  } while (++iter_count < max_iter);
  
  // Optimized final coefficient computation
  arma::uvec active_predictors = arma::find(new_betas != 0);
  arma::uvec active_samples = arma::find(new_trim == 0);
  
  if (!active_predictors.empty() && !active_samples.empty()) {
    arma::mat X_active = x_subset.submat(active_samples, active_predictors);
    arma::vec y_active = y_sc.elem(active_samples);
    
    arma::mat XtX_active = X_active.t() * X_active;
    arma::vec Xty_active = X_active.t() * y_active;
    XtX_active.diag() += 1e-8;
    
    arma::vec beta_active = arma::solve(XtX_active, Xty_active, arma::solve_opts::fast);
    new_betas.elem(active_predictors) = beta_active;
  }
  
  coef_mat_candidate.col(group).zeros();
  group_vec(0) = group;
  coef_mat_candidate.submat(model_subspace, group_vec) = new_betas;
  
  // Updating subset indices for the group
  Update_Subset_Indices_Candidate(group);
  
  // Updating the active samples for the group
  Update_Active_Samples_Candidate(group, new_trim);
}

// OPTIMIZED: Projection functions using partial sorting
void EnsembleModel::Project_Coef(arma::vec& coef_vector) {
  
  if (t >= coef_vector.n_elem) return;  // No projection needed
  
  // Use nth_element for O(n) partial sorting instead of full O(n log n) sort
  arma::vec abs_coef = arma::abs(coef_vector);
  std::vector<double> abs_vec = arma::conv_to<std::vector<double>>::from(abs_coef);
  
  if (t < abs_vec.size()) {
    std::nth_element(abs_vec.begin(), abs_vec.begin() + t, abs_vec.end(), std::greater<double>());
    double threshold = abs_vec[t];
    
    // Zero out elements below threshold
    coef_vector.elem(arma::find(abs_coef < threshold)).zeros();
  }
}

void EnsembleModel::Project_Trim(arma::vec& trim_vector) {
  
  if (h >= n) return;  // No trimming needed
  
  arma::vec abs_trim = arma::abs(trim_vector);
  std::vector<double> abs_vec = arma::conv_to<std::vector<double>>::from(abs_trim);
  
  if ((n-h) < abs_vec.size()) {
    std::nth_element(abs_vec.begin(), abs_vec.begin() + (n-h), abs_vec.end(), std::greater<double>());
    double threshold = abs_vec[n-h];
    
    trim_vector.elem(arma::find(abs_trim < threshold)).zeros();
  }
}

double EnsembleModel::Compute_Group_Loss(arma::mat& x, arma::vec& y, arma::vec& betas, arma::vec& trim) {
  
  return arma::sum(arma::square(y - x * betas - trim));
}

void EnsembleModel::Update_Subset_Indices(arma::uword& group) {
  
  subset_active.zeros();
  subset_active(arma::find(coef_mat.col(group) != 0)).ones();
  subset_indices.col(group) = subset_active;
  invalidate_cache();  // Row sums changed
}
void EnsembleModel::Update_Subset_Indices_Candidate(arma::uword& group) {
  
  subset_active.zeros();
  subset_active(arma::find(coef_mat_candidate.col(group) != 0)).ones();
  subset_indices_candidate.col(group) = subset_active;
  invalidate_candidate_cache();
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
  
  arma::mat predictions = x * final_coef;  
  for(arma::uword group = 0; group < n_models; group++)
    models_loss(group) = arma::mean(arma::square(y - final_intercept(group) - predictions.col(group)));
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
    
    // Update main cache when accepting candidate
    invalidate_cache();
  }
}

// OPTIMIZED: Subspace computation with caching
arma::uvec EnsembleModel::Get_Model_Subspace(arma::uword& group) {
  update_row_sums_cache();
  return arma::find((row_sums_subset - subset_indices.col(group)) < u);
}
arma::uvec EnsembleModel::Get_Model_Subspace_Candidate(arma::uword& group) {
  update_row_sums_candidate_cache();
  return arma::find((row_sums_subset_candidate - subset_indices_candidate.col(group)) < u);
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
  
  return (y_test - arma::mean(final_intercept) - arma::mean(x_test * final_coef, 1));
}
arma::vec EnsembleModel::Prediction_Square_Residuals_Ensemble(arma::mat& x_test, arma::vec& y_test) {
  
  return arma::square(y_test - arma::mean(final_intercept) - arma::mean(x_test * final_coef, 1));
}
arma::vec EnsembleModel::Prediction_Square_Residuals_Models(arma::mat& x_test, arma::vec& y_test) {
  
  arma::vec prediction_residuals = arma::zeros(y_test.n_elem);
  for (arma::uword group = 0; group < n_models; group++) {
    
    prediction_residuals += arma::square(y_test - final_intercept(group) - x_test * final_coef.col(group));
  }
  return prediction_residuals/n_models;
}