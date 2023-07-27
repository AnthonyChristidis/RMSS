#' 
#' @useDynLib RMSS
#' @importFrom Rcpp sourceCpp
#' 
#' @importFrom stats coef median mad 
#'
#' @title Robust Multi-Model Subset Selection
#' 
#' @description \code{RMSS} performs robust multi-model subset selection.
#' 
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_models Number of models into which the variables are split.
#' @param h_grid Grid for robustness parameter.
#' @param t_grid Grid for sparsity parameter.
#' @param u_grid Grid for diversity parameter.
#' @param initial_estimator Method used for initial estimator. Must be one of "srlars" (default) or "robStepSplitReg".
#' @param tolerance Tolerance level for convergence of PSBGD algorithm.
#' @param max_iter Maximum number of iterations in PSBGD algorithm.
#' @param neighborhood_search Neighborhood search to improve solution. Default is FALSE.
#' @param neighborhood_search_tolerance Tolerance parameter for neighborhood search. Default is 1e-1.
#' 
#' @return An object of class RMSS
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{coef.RMSS}}, \code{\link{predict.RMSS}}
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
#' # RMSS
#' rmss_fit <- RMSS(x = x_train, y = y_train,
#'                  n_models = 3,
#'                  h_grid = c(35), t_grid = c(6, 8, 10), u_grid = c(1:3),
#'                  initial_estimator = "srlars",
#'                  tolerance = 1e-1,
#'                  max_iter = 1e3,
#'                  neighborhood_search = FALSE,
#'                  neighborhood_search_tolerance = 1e-1)
#' rmss_coefs <- coef(rmss_fit, 
#'                    h_ind = 1, t_ind = 2, u_ind = 1,
#'                    group_index = 1:rmss_fit$n_models)
#' sens_rmss <- sum(which((rmss_coefs[-1]!=0)) <= p.active)/p.active
#' spec_rmss <- sum(which((rmss_coefs[-1]!=0)) <= p.active)/sum(rmss_coefs[-1]!=0)
#' rmss_preds <- predict(rmss_fit, newx = x_test,
#'                       h_ind = 1, t_ind = 2, u_ind = 1,
#'                       group_index = 1:rmss_fit$n_models,
#'                       dynamic = FALSE)
#' rmss_mspe <- mean((y_test - rmss_preds)^2)/sigma^2
#' 
RMSS <- function(x, y,
                 n_models,
                 h_grid, t_grid, u_grid,
                 initial_estimator = c("srlars", "robStepSplitReg")[1],
                 tolerance = 1e-1,
                 max_iter = 1e3,
                 neighborhood_search = FALSE,
                 neighborhood_search_tolerance = 1e-1){
  
  # Check input data
  DataCheck(x, y,
            n_models,
            h_grid, t_grid, u_grid,
            initial_estimator,
            tolerance,
            max_iter,
            neighborhood_search,
            neighborhood_search_tolerance)
  
  # Shuffle the data
  n <- nrow(x)
  p <- ncol(x)
  random.permutation <- sample(1:n, n)
  x <- x[random.permutation, ]
  y <- y[random.permutation]
  
  # Initial split of predictors
  if(initial_estimator == "srlars"){
    
    initial_selections <- srlars::srlars(x, y,
                                         n_models = n_models,
                                         model_saturation = "fixed",
                                         model_size = min(n - 1, floor(p/n_models)),
                                         robust = TRUE,
                                         compute_coef = FALSE)$selections
    
  } else if(initial_estimator == "robStepSplitReg"){
    
    initial_selections <- robStepSplitReg::robStepSplitReg(x, y, 
                                                           n_models = n_models,
                                                           model_saturation = "fixed",
                                                           model_size = min(n - 1, floor(p/n_models)),
                                                           robust = TRUE,
                                                           compute_coef = FALSE)$selections
  }
  initial_split <- matrix(0, nrow = p, ncol = n_models)
  for(model_id in 1:n_models)
    initial_split[initial_selections[[model_id]], model_id] <- 1
  
  # CPP parameter for neighborhood search
  neighborhood_search_cpp <- as.numeric(neighborhood_search)
  
  # Invoking the CPP code for RMSS
  output <- RInterface(x, y,
                       n_models,
                       h_grid, t_grid, u_grid,
                       tolerance,
                       max_iter,
                       initial_split,
                       neighborhood_search_cpp,
                       neighborhood_search_tolerance)
  
  # Formatting of coefficients
  for(h_ind in 1:length(h_grid))
    for(t_ind in 1:length(t_grid))
      for(u_ind in 1:length(u_grid)){
        
        output$active_samples[[h_ind]][[t_ind]][[u_ind]] <- matrix(output$active_samples[[h_ind]][[t_ind]][[u_ind]], ncol = n_models, byrow = FALSE)[order(random.permutation),]
        output$coef[[h_ind]][[t_ind]][[u_ind]] <- matrix(output$coef[[h_ind]][[t_ind]][[u_ind]], ncol = n_models, byrow = FALSE)
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
  class(output) <- append("RMSS", class(output))
  
  # Returning the output from the stepwise algorithm
  return(output)
}



