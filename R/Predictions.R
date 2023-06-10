#' 
#' @title Predictions for RMSS Object
#' 
#' @description \code{predict.RMSS} returns the predictions for a RMSS object.
#' 
#' @method predict RMSS
#' 
#' @param object An object of class RMSS.
#' @param newx New data for predictions.
#' @param h_ind Index for robustness parameter.
#' @param t_ind Index for sparsity parameter.
#' @param u_ind Index for diversity parameter.
#' @param group_index Groups included in the ensemble. Default setting includes all the groups.
#' @param dynamic Argument to determine whether dynamic predictions are used based on deviating cells. Default is FALSE.
#' @param ... Additional arguments for compatibility.
#' 
#' @return The predictions for the RMSS object.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{RMSS}}
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
predict.RMSS <- function(object, newx, 
                         h_ind, t_ind, u_ind,
                         group_index = NULL, 
                         dynamic = FALSE,
                         ...){
  
  if(!dynamic){
    
    ensemble.coef <- coef(object, h_ind, t_ind, u_ind, group_index = group_index)
    output <- ensemble.coef[1] + as.numeric(newx %*% ensemble.coef[-1])
    return(output)
  } else{
    
    DDC_cells <- cellWise::DDCpredict(newx, object$DDCx)$indcells
    x_test_cells <- newx
    x_test_cells[DDC_cells] <- NA
    cells_id <- apply(x_test_cells, 1, function(x) return(which(is.na(x))))
    var_selections <- apply(object$coef[[h_ind]][[t_ind]][[u_ind]], 2,
                            function(x) return(which(x!=0)))
    selected_models <- matrix(0, nrow = nrow(newx), ncol = object$n_models)
    for(model_id in 1:object$n_models){
      
      selected_models[, model_id] <- sapply(cells_id, function(x) return(any(x %in% var_selections[, model_id])), 
                                            simplify = TRUE)
    }
    selected_models <- lapply(1:nrow(newx), function(x, selected_models) return(which(selected_models[x, ] != 0)),
                              selected_models = selected_models)
    output <- numeric(nrow(newx))
    for(obs_id in 1:nrow(newx)){
      
      ensemble.coef <- coef(object, 
                            h_ind = h_ind, t_ind = t_ind, u_ind = u_ind,
                            group_index = selected_models[[obs_id]])
      output[obs_id] <- ensemble.coef[1] + as.numeric(newx[obs_id,] %*% ensemble.coef[-1])
    }
    return(output)
  }
}
#' 
#' @title Predictions for cv.RMSS Object
#' 
#' @description \code{predict.cv.RMSS} returns the predictions for a cv.RMSS object.
#' 
#' @method predict cv.RMSS
#' 
#' @param object An object of class cv.RMSS.
#' @param newx New data for predictions.
#' @param h_ind Index for robustness parameter.
#' @param t_ind Index for sparsity parameter.
#' @param u_ind Index for diversity parameter.
#' @param group_index Groups included in the ensemble. Default setting includes all the groups.
#' @param dynamic Argument to determine whether dynamic predictions are used based on deviating cells. Default is FALSE.
#' @param ... Additional arguments for compatibility.
#' 
#' @return The predictions for the cv.RMSS object.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{cv.RMSS}}
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
predict.cv.RMSS <- function(object, newx, 
                            h_ind = NULL, t_ind = NULL, u_ind = NULL,
                            group_index = NULL, 
                            dynamic = FALSE,
                            ...){
  
  if(is.null(h_ind))
    h_ind <- object$h_opt
  if(is.null(t_ind))
    t_ind <- object$t_opt
  if(is.null(u_ind))
    u_ind <- object$u_opt
  
  if(!dynamic){
    
    ensemble.coef <- coef(object, h_ind, t_ind, u_ind, group_index = group_index)
    output <- ensemble.coef[1] + as.numeric(newx %*% ensemble.coef[-1])
    return(output)
  } else{
    
    DDC_cells <- cellWise::DDCpredict(newx, object$DDCx)$indcells
    x_test_cells <- newx
    x_test_cells[DDC_cells] <- NA
    cells_id <- apply(x_test_cells, 1, function(x) return(which(is.na(x))))
    var_selections <- apply(object$coef[[h_ind]][[t_ind]][[u_ind]], 2,
                            function(x) return(which(x!=0)))
    selected_models <- matrix(0, nrow = nrow(newx), ncol = object$n_models)
    for(model_id in 1:object$n_models){
      
      selected_models[, model_id] <- sapply(cells_id, function(x) return(any(x %in% var_selections[, model_id])), 
                                            simplify = TRUE)
    }
    selected_models <- lapply(1:nrow(newx), function(x, selected_models) return(which(selected_models[x, ] != 0)),
                              selected_models = selected_models)
    empty_models <- which(lapply(selected_models, length) == 0)
    for(empty_id in empty_models)
      selected_models[[empty_id]] <- 1:object$n_models
    output <- numeric(nrow(newx))
    for(obs_id in 1:nrow(newx)){
      
      ensemble.coef <- coef(object, 
                            h_ind = h_ind, t_ind = t_ind, u_ind = u_ind,
                            group_index = selected_models[[obs_id]])
      output[obs_id] <- ensemble.coef[1] + as.numeric(newx[obs_id,] %*% ensemble.coef[-1])
    }
    return(output)
  }
}


