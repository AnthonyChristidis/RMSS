# ---------------------------------
# Function to replicate RCPP folds
# ---------------------------------

ReplicateRCPPFolds <- function(n, n_folds){
  
  fold_int <- seq(from = 1, to = n, by = floor(n/n_folds)) 
  if(floor(n/n_folds)*n_folds == n)
    fold_int <- c(fold_int, n + 1) else{
      
      excess <- n - (fold_int[length(fold_int)] - 1)
      for(adjust in (n_folds + 1 - (excess - 1)):(n_folds + 1))
        fold_int[adjust:(n_folds + 1)] <- fold_int[adjust:(n_folds + 1)] + 1
    }
  fold_list <- list()
  for(fold in 1:n_folds)
    fold_list[[fold]] <- fold_int[fold]:(fold_int[fold + 1] - 1)
  fold_final <- list()
  for(fold in 1:n_folds)
    fold_final[[fold]] <- unlist(fold_list)[-fold_list[[fold]]]
  return(fold_final)
}


