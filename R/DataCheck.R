# --------------------------------------
# Checking Input Data for RMSS Function
# --------------------------------------
DataCheck <- function(x, y,
                      n_models,
                      h_grid, t_grid, u_grid,
                      initial_estimator,
                      tolerance,
                      max_iter,
                      neighborhood_search,
                      neighborhood_search_tolerance){
  
  # Checking the input for the design matrix (x) and the response vector (y)
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame.")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric.")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values.")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values.")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector.")
      }
      # Force to vector if input was a matrix
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows.")
    }
  }
  
  # Checking the input for the number of models
  if(!is.null(n_models)){
    if (!inherits(n_models, "numeric")) {
      stop("n_models should be numeric.")
    } else if (any(!n_models == floor(n_models), n_models <= 0, n_models > nrow(x))) {
      stop("n_models should be a positive integer greater than 0 and less than the sample size.")
    }
  }
  
  # Checking the input for the trimming grid
  if(!is.null(h_grid)){
    if (!inherits(h_grid, c("numeric", "integer"))) {
      stop("h_grid should be numeric or an integer vector.")
    } else if (any(!(h_grid == floor(h_grid)), h_grid < 1, h_grid > nrow(x))) {
      stop("h_grid should be a positive integer no larger than the sample size.")
    }
  }
  
  # Checking the input for the sparsity grid
  if(!is.null(t_grid)){
    if (!inherits(t_grid, c("numeric", "integer"))) {
      stop("t_grid should be numeric or an integer vector.")
    } else if (any(!(t_grid == floor(t_grid)), t_grid < 1, t_grid >= nrow(x))) {
      stop("t_grid should be a positive integer less than the sample size.")
    }
  }
  
  # Checking the input for the diversity grid
  if(!is.null(u_grid)){
    if (!inherits(u_grid, c("numeric", "integer"))) {
      stop("u_grid should be numeric or an integer vector.")
    } else if (any(!(u_grid == floor(u_grid)), u_grid < 1, u_grid > n_models)) {
      stop("u_grid should be a positive integer less than the number of models.")
    }
  }
  
  # Check input for initial estimator 
  if(!(initial_estimator %in% c("srlars", "robStepSplitReg")))
    stop("initial_estimator should be one of \"srlars\" or \"robStepSplitReg\".")
  
  # Checking input for the tolerance parameter
  if(!is.null(tolerance))
    if (!inherits(tolerance, "numeric")) {
      stop("tolerance should be numeric.")
    } else if (length(tolerance) != 1) {
      stop("tolerance should be a single numeric value.")
    }
  
  # Checking the input for the maximum number of iterations
  if(!is.null(max_iter)){
    if (!inherits(max_iter, "numeric")) {
      stop("max_iter should be numeric.")
    } else if (any(!max_iter == floor(max_iter), max_iter <= 0)) {
      stop("max_iter should be a positive integer greater than 0.")
    }
  }
  
  # Check neighborhood search parameter
  if(!(neighborhood_search %in% c(TRUE, FALSE)))
    stop("neighborhood_search should be TRUE or FALSE.")
  
  # Checking input for the neighborhood_search_tolerance parameter
  if(!is.null(neighborhood_search_tolerance))
    if (!inherits(neighborhood_search_tolerance, "numeric")) {
      stop("neighborhood_search_tolerance should be numeric.")
    } else if (length(neighborhood_search_tolerance) != 1) {
      stop("neighborhood_search_tolerance should be a single numeric value.")
    }
}

# -----------------------------------------
# Checking Input Data for cv.RMSS Function
# -----------------------------------------
DataCheckCV <- function(x, y,
                        n_models,
                        h_grid, t_grid, u_grid,
                        initial_estimator,
                        tolerance,
                        max_iter,
                        neighborhood_search,
                        neighborhood_search_tolerance,
                        n_folds,
                        cv_criterion,
                        alpha,
                        gamma, 
                        n_threads){
  
  # Checking the input for the design matrix (x) and the response vector (y)
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame.")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric.")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values.")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values.")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector.")
      }
      # Force to vector if input was a matrix
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows.")
    }
  }
  
  # Checking the input for the number of models
  if(!is.null(n_models)){
    if (!inherits(n_models, "numeric")) {
      stop("n_models should be numeric.")
    } else if (any(!n_models == floor(n_models), n_models <= 0, n_models > nrow(x))) {
      stop("n_models should be a positive integer greater than 0 and less than the sample size.")
    }
  }
  
  # Checking the input for the trimming grid
  if(!is.null(h_grid)){
    if (!inherits(h_grid, c("numeric", "integer"))) {
      stop("h_grid should be numeric or an integer vector.")
    } else if (any(!(h_grid == floor(h_grid)), h_grid < 1, h_grid > nrow(x))) {
      stop("h_grid should be a positive integer no larger than the sample size.")
    }
  }
  
  # Checking the input for the sparsity grid
  if(!is.null(t_grid)){
    if (!inherits(t_grid, c("numeric", "integer"))) {
      stop("t_grid should be numeric or an integer vector.")
    } else if (any(!(t_grid == floor(t_grid)), t_grid < 1, t_grid >= nrow(x))) {
      stop("t_grid should be a positive integer less than the sample size.")
    }
  }
  
  # Checking the input for the diversity grid
  if(!is.null(u_grid)){
    if (!inherits(u_grid, c("numeric", "integer"))) {
      stop("u_grid should be numeric or an integer vector.")
    } else if (any(!(u_grid == floor(u_grid)), u_grid < 1, u_grid > n_models)) {
      stop("u_grid should be a positive integer less than the number of models.")
    }
  }
  
  # Check input for initial estimator 
  if(!(initial_estimator %in% c("srlars", "robStepSplitReg")))
    stop("initial_estimator should be one of \"srlars\" or \"robStepSplitReg\".")
  
  # Checking input for the tolerance parameter
  if(!is.null(tolerance))
    if (!inherits(tolerance, "numeric")) {
      stop("tolerance should be numeric.")
    } else if (length(tolerance) != 1) {
      stop("tolerance should be a single numeric value.")
    }
  
  # Checking the input for the maximum number of iterations
  if(!is.null(max_iter)){
    if (!inherits(max_iter, "numeric")) {
      stop("max_iter should be numeric.")
    } else if (any(!max_iter == floor(max_iter), max_iter <= 0)) {
      stop("max_iter should be a positive integer greater than 0.")
    }
  }
  
  # Check neighborhood search parameter
  if(!(neighborhood_search %in% c(TRUE, FALSE)))
    stop("neighborhood_search should be TRUE or FALSE.")
  
  # Checking input for the neighborhood_search_tolerance parameter
  if(!is.null(neighborhood_search_tolerance))
    if (!inherits(neighborhood_search_tolerance, "numeric")) {
      stop("neighborhood_search_tolerance should be numeric.")
    } else if (length(neighborhood_search_tolerance) != 1) {
      stop("neighborhood_search_tolerance should be a single numeric value.")
    }
  
  # Check input for number of folds
  if(!inherits(n_folds, "numeric")) {
    stop("n_folds should be numeric")
  } else if(any(!n_folds == floor(n_folds), n_folds <= 0)) {
    stop("n_folds should be a positive integer")
  }
  
  # Check CV criterion parameter
  if(!(cv_criterion %in% c("tau", "trimmed")))
    stop("cv_criterion should be one of \"tau\" or \"trimmed\".")
  
  # Check alpha value
  if(!inherits(alpha, "numeric")) {
    stop("alpha should be numeric.")
  } else if(!all(alpha <= 1, alpha > 0)) {
    stop("alpha should be between 0 and 1.")
  }
  
  # Check gamma value
  if(!inherits(gamma, "numeric")) {
    stop("gamma should be numeric.")
  } else if(!all(gamma <= 1, gamma > 0)) {
    stop("gamma should be between 0 and 1.")
  }
  
  # Check input for number of threads
  if(!inherits(n_threads, "numeric")) {
    stop("n_threads should be numeric")
  } else if(any(!n_threads == floor(n_threads), n_threads <= 0)) {
    stop("n_threads should be a positive integer")
  }
}

# ---------------------------------------
# Checking Input Data for coef Functions
# ---------------------------------------
DataCheckCoef <- function(object,
                          h_ind, t_ind, u_ind,
                          individual_models){
  
  # Check index of trimming parameter
  if(!is.null(h_ind)){
    if(any(!inherits(h_ind, c("numeric", "integer")), length(h_ind) > 1))
      stop("h_ind should be a single numeric or integer value.") else if(!(h_ind %in% c(1:length(object$h_grid))))
        stop("h_ind is out of range.")
  }
  
  # Check index of sparsity parameter
  if(!is.null(t_ind)){
    if(any(!inherits(t_ind, c("numeric", "integer")), length(t_ind) > 1))
      stop("t_ind should be a single numeric or integer value.") else if(!(t_ind %in% c(1:length(object$t_grid))))
        stop("t_ind is out of range.")
  }
  
  # Check index of diversity parameter
  if(!is.null(u_ind)){
    if(any(!inherits(u_ind, c("numeric", "integer")), length(u_ind) > 1))
      stop("u_ind should be a single numeric or integer value.") else if(!(u_ind %in% c(1:length(object$u_grid))))
        stop("u_ind is out of range.")
  }
  
  # Check individual models parameter
  if(!(individual_models %in% c(TRUE, FALSE)))
    stop("individual_models should be TRUE or FALSE.")
}

# ------------------------------------------
# Checking Input Data for predict Functions
# ------------------------------------------
DataCheckPredict <- function(object, newx,
                             h_ind, t_ind, u_ind){
  
  # Check number of predictors in test data
  if(inherits(newx, c("matrix", "data.frame"))){
    if(ncol(newx) != object$p)
      stop("The number of predictors in the new data does not match the number of predictors in the training data.")
  }
  
  # Check index of trimming parameter
  if(!is.null(h_ind)){
    if(any(!inherits(h_ind, c("numeric", "integer")), length(h_ind) > 1))
      stop("h_ind should be a single numeric or integer value.") else if(!(h_ind %in% c(1:length(object$h_grid))))
        stop("h_ind is out of range.")
  }
  
  # Check index of sparsity parameter
  if(!is.null(t_ind)){
    if(any(!inherits(t_ind, c("numeric", "integer")), length(t_ind) > 1))
      stop("t_ind should be a single numeric or integer value.") else if(!(t_ind %in% c(1:length(object$t_grid))))
        stop("t_ind is out of range.")
  }
  
  # Check index of diversity parameter
  if(!is.null(u_ind)){
    if(any(!inherits(u_ind, c("numeric", "integer")), length(u_ind) > 1))
      stop("u_ind should be a single numeric or integer value.") else if(!(u_ind %in% c(1:length(object$u_grid))))
        stop("u_ind is out of range.")
  }
}

