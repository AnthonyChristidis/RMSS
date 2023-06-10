#' 
#' @title Cross-Validatoin for Robust Multi-Model Subset Selection
#' 
#' @description \code{cv.RMSS} performs the cross-validation procedure for robust multi-model subset selection.
#' 
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_models Number of models into which the variables are split.
#' @param h_grid Grid for robustness parameter.
#' @param t_grid Grid for sparsity parameter.
#' @param u_grid Grid for diversity parameter.
#' @param tolerance Tolerance level for convergence of PSBGD algorithm.
#' @param max_iter Maximum number of iterations in PSBGD algorithm.
#' @param neighborhood_search Neighborhood search to improve solution. Default is FALSE.
#' @param neighborhood_search_tolerance Tolerance parameter for neighborhood search. Default is 1e-1.
#' @param n_folds Number of folds for cross-validation procedure. Default is 5.
#' @param alpha Proportion of trimmed samples for cross-validation procedure. Default is 1/4.
#' @param gamma Weight parameter for ensemble MSPE (gamma) and average MSPE of individual models (1 - gamma). Default is 1.
#' @param n_threads Number of threads used by OpenMP for multithreading over the folds. Default is 1.
#'  
#' @return An object of class cv.RMSS
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{coef.cv.RMSS}}, \code{\link{predict.cv.RMSS}}
#' 
#' @examples
#' # Simulation parameters
#' n <- 50
#' p <- 100
#' rho <- 0.8
#' rho.inactive <- 0.2
#' group.size <- 5
#' p.active <- 15
#' snr <- 2
#' contamination.prop <- 0.3
#' 
#' # Setting the seed
#' set.seed(0)
#' 
#' # Block Correlation
#' sigma.mat <- matrix(0, p, p)
#' sigma.mat[1:p.active, 1:p.active] <- rho.inactive
#' for(group in 0:(p.active/group.size - 1))
#'   sigma.mat[(group*group.size+1):(group*group.size+group.size),
#'             (group*group.size+1):(group*group.size+group.size)] <- rho
#' diag(sigma.mat) <- 1
#' 
#' # Simulation of beta vector
#' true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), 
#'                rep(0, p - p.active))
#' 
#' # Setting the SD of the variance
#' sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))
#' 
#' # Simulation of test data
#' m <- 2e3
#' x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
#' y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)
#' 
#' # Simulation of uncontaminated data 
#' x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
#' y <- x %*% true.beta + rnorm(n, 0, sigma)
#' 
#' # Contamination of data 
#' contamination_indices <- 1:floor(n*contamination.prop)
#' k_lev <- 2
#' k_slo <- 100
#' x_train <- x
#' y_train <- y
#' beta_cont <- true.beta
#' beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
#' beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
#' for(cont_id in contamination_indices){
#'  
#'  a <- runif(p, min = -1, max = 1)
#'  a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
#'   x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + k_lev * a / 
#'                         as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
#'   y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
#' }
#' 
#' # CV RMSS
#' rmss_fit <- cv.RMSS(x = x_train, y = y_train,
#'                     n_models = 3,
#'                     h_grid = c(35), t_grid = c(6, 8, 10), u_grid = c(1:3),
#'                     tolerance = 1e-1,
#'                     max_iter = 1e3,
#'                     neighborhood_search = FALSE,
#'                     neighborhood_search_tolerance = 1e-1,
#'                     n_folds = 5,
#'                     alpha = 1/4,
#'                     gamma = 1, 
#'                     n_threads = 1)
#' rmss_coefs <- coef(rmss_fit, 
#'                    h_ind = rmss_fit$h_opt, 
#'                    t_ind = rmss_fit$t_opt, 
#'                    u_ind = rmss_fit$u_opt,
#'                    group_index = 1:rmss_fit$n_models)
#' sens_rmss <- sum(which((rmss_coefs[-1]!=0)) <= p.active)/p.active
#' spec_rmss <- sum(which((rmss_coefs[-1]!=0)) <= p.active)/sum(rmss_coefs[-1]!=0)
#' rmss_preds <- predict(rmss_fit, newx = x_test,
#'                       h_ind = rmss_fit$h_opt, 
#'                       t_ind = rmss_fit$t_opt, 
#'                       u_ind = rmss_fit$u_opt,
#'                       group_index = 1:rmss_fit$n_models,
#'                       dynamic = FALSE)
#' rmss_mspe <- mean((y_test - rmss_preds)^2)/sigma^2
#'  
cv.RMSS <- function(x, y,
                    n_models,
                    h_grid, t_grid, u_grid,
                    tolerance = 1e-1,
                    max_iter = 1e3,
                    neighborhood_search = FALSE,
                    neighborhood_search_tolerance = 1e-1,
                    n_folds = 5,
                    alpha = 1/4,
                    gamma = 1, 
                    n_threads = 1){
  
  # Shuffle the data
  n <- nrow(x)
  p <- ncol(x)
  random.permutation <- sample(1:n, n)
  x <- x[random.permutation, ]
  y <- y[random.permutation]
  
  # Initial split of predictors
  initial_selections <- robStepSplitReg::robStepSplitReg(x, y, 
                                                         n_models = n_models,
                                                         model_saturation = "fixed",
                                                         model_size = n - 1,
                                                         robust = TRUE,
                                                         compute_coef = FALSE)$selections
  initial_split <- matrix(0, nrow = p, ncol = n_models)
  for(model_id in 1:n_models)
    initial_split[initial_selections[[model_id]], model_id] <- 1
  
  # Creation of the folds (CPP replication)
  cpp_folds <- ReplicateRCPPFolds(n, n_folds)
  
  # Array for initial selections
  initial_split_array <- array(dim = c(nrow(initial_split), ncol(initial_split), n_folds))

  # Filling arrays for data scaling and initial selections
  for(fold in 1:n_folds){
    
    # Initial split of predictors
    initial_selections_fold <- robStepSplitReg::robStepSplitReg(x[cpp_folds[[fold]],], y[cpp_folds[[fold]]], 
                                                                n_models = n_models,
                                                                model_saturation = "fixed",
                                                                model_size = length(cpp_folds[[fold]]) - 1,
                                                                robust = TRUE,
                                                                compute_coef = FALSE)$selections
    initial_split_array[,, fold] <- matrix(0, nrow = p, ncol = n_models)
    for(model_id in 1:n_models)
      initial_split_array[initial_selections_fold[[model_id]], model_id, fold] <- 1
  }
  
  # CPP parameter for neighborhood search
  neighborhood_search_cpp <- as.numeric(neighborhood_search)
  
  # Invoking the CPP code for RMSS
  output <- RInterfaceCV(x, y,
                         n_models,
                         h_grid, t_grid, u_grid,
                         tolerance,
                         max_iter,
                         initial_split_array,
                         initial_split,
                         neighborhood_search_cpp,
                         neighborhood_search_tolerance,
                         n_folds,
                         alpha,
                         gamma,
                         n_threads)
  
  # Formatting the list of output
  for(h_ind in 1:length(h_grid)){
    for(t_ind in 1:length(t_grid)){
      
      output$cv_error[[h_ind]][[t_ind]] <- as.list(output$cv_error[[h_ind]][[t_ind]])
      
      for(u_ind in 1:length(u_grid)){
        
        output$active_samples[[h_ind]][[t_ind]][[u_ind]] <- matrix(output$active_samples[[h_ind]][[t_ind]][[u_ind]], ncol = n_models, byrow = FALSE)[order(random.permutation),]
        output$coef[[h_ind]][[t_ind]][[u_ind]] <- matrix(output$coef[[h_ind]][[t_ind]][[u_ind]], ncol = n_models, byrow = FALSE)
      }
    }
  }
  
  # Formatting CV error
  cv_error <- lapply(1:n_models, function(t) return(matrix(NA, nrow = length(h_grid), ncol = length(t_grid))))
  for(h_ind in 1:length(h_grid))
    for(t_ind in 1:length(t_grid))
      for(u_ind in 1:length(u_grid))
        cv_error[[u_ind]][h_ind, t_ind] <- output$cv_error[[h_ind]][[t_ind]][[u_ind]]
  output$cv_error <- cv_error
  
  # Optimal parameters (RMSS)
  output$u_opt <- which.min(sapply(output$cv_error, min))
  output$h_opt <- which(output$cv_error[[output$u_opt]] == min(output$cv_error[[output$u_opt]]), arr.ind = TRUE)[1]
  output$t_opt <- which(output$cv_error[[output$u_opt]] == min(output$cv_error[[output$u_opt]]), arr.ind = TRUE)[2]
  
  # Optimal parameters (RBSS)
  if(u_grid[length(u_grid)] == n_models){
    
    output$rbss_h_opt <- which(output$cv_error[[n_models]] == min(output$cv_error[[n_models]]), arr.ind = TRUE)[1]
    output$rbss_t_opt <- which(output$cv_error[[n_models]] == min(output$cv_error[[n_models]]), arr.ind = TRUE)[2]
  }
  
  # Adding specifications of output
  output$n_models <- n_models
  output$h_grid <- h_grid
  output$t_grid <- t_grid
  output$u_grid <- u_grid
  output$tolerance <- tolerance
  output$max_iter <- max_iter
  output$n <- nrow(x)
  output$p <- ncol(x)
  output$DDCx <- cellWise::DDC(x, DDCpars = list(fastDDC = TRUE, nbngbrs = p-1, silent = TRUE))
  
  # Create the object of class "stepSplitReg"
  class(output) <- append("cv.RMSS", class(output))
  
  # Returning the output from the stepwise algorithm
  return(output)
}



