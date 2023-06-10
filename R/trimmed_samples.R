#' 
#' @title Trimmed samples for RMSS or cv.RMSS Object
#' 
#' @description \code{trimmed_samples} returns the coefficients for a RMSS or cv.RMSS object.
#' 
#' @param object An object of class RMSS
#' @param h_ind Index for robustness parameter.
#' @param t_ind Index for sparsity parameter.
#' @param u_ind Index for diversity parameter.
#' @param group_index Groups included in the ensemble. Default setting includes all the groups.
#' @param ... Additional arguments for compatibility.
#' 
#' @return The trimmed samples for the RMSS or cv.RMSS object.
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
#' trimmed_id <- trimmed_samples(rmss_fit, h_ind = 1, t_ind = 1, u_ind = 1)
#' 
trimmed_samples <- function(object, 
                            h_ind = NULL, t_ind = NULL, u_ind = NULL,
                            group_index = NULL, 
                            ...){
  
  if(!(any(class(object) %in% c("RMSS", "cv.RMSS"))))
    stop("The object is not of class \"RMSS\" or \"cv.RMSS\".")
  
  if(any(class(object) == "cv.RMSS")){
    
    if(is.null(h_ind))
      h_ind <- object$h_opt
    if(is.null(t_ind))
      t_ind <- object$t_opt
    if(is.null(u_ind))
      u_ind <- object$u_opt
    
  } else if(any(class(object) == "RMSS")){
    
    if(any(is.null(h_ind), is.null(t_ind), is.null(u_ind)))
      stop("The arguments \"h_ind\", \"t_ind\" and \"u_ind\" must be specified.")
  }
  
  if(is.null(group_index)){
    
    return(apply(object$active_samples[[h_ind]][[t_ind]][[u_ind]], 2, function(x) return(which(x == 0))))
    
  } else{
    
    if(any(!(group_index %in% 1:object$n_models)))
      stop("The group index is invalid.")
    
    return(apply(object$active_samples[[h_ind]][[t_ind]][[u_ind]], 2, function(x) return(which(x == 0)))[, group_index])
  }
}

