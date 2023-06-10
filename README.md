[![CRAN checks](https://badges.cranchecks.info/summary/RMSS.svg)](https://cran.r-project.org/web/checks/check_results_RMSS.html) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/RMSS)](https://cran.r-project.org/package=RMSS) [![Downloads](http://cranlogs.r-pkg.org/badges/RMSS)](https://cran.r-project.org/package=RMSS)

RMSS
=====

This package provides functions for performing robust multi-model subset selection.

------------------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=RMSS).

```{r installation, eval = FALSE}
install.packages("RMSS", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/RMSS)

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/RMSS")
```

### Usage

``` r
# Required libraries
install.packages("mvnfast")

# Simulation parameters
n <- 50
p <- 500
rho <- 0.8
rho.inactive <- 0.2
group.size <- 25
p.active <- 100
snr <- 1
contamination.prop <- 0.2

# Setting the seed
set.seed(0)

# Block Correlation
sigma.mat <- matrix(0, p, p)
sigma.mat[1:p.active, 1:p.active] <- rho.inactive
for(group in 0:(p.active/group.size - 1))
  sigma.mat[(group*group.size+1):(group*group.size+group.size),(group*group.size+1):(group*group.size+group.size)] <- rho
diag(sigma.mat) <- 1

# Simulation of beta vector
true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))

# Setting the SD of the variance
sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))

# Simulation of test data
m <- 2e3
x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)

# Simulation of uncontaminated data 
x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
y <- x %*% true.beta + rnorm(n, 0, sigma)

# Contamination of data 
contamination_indices <- 1:floor(n*contamination.prop)
k_lev <- 2
k_slo <- 100
x_train <- x
y_train <- y
beta_cont <- true.beta
beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
for(cont_id in contamination_indices){
  
  a <- runif(p, min = -1, max = 1)
  a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
  x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
  y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
}

# CV RMSS
cv.rmss_fit <- cv.RMSS(x = x_train, y = y_train,
                       n_models = 10,
                       h_grid = c(35, 40), t_grid = c(8, 10, 12), u_grid = c(1:10),
                       tolerance = 1e-1,
                       max_iter = 1e3,
                       neighborhood_search = FALSE,
                       neighborhood_search_tolerance = 1e-1,
                       n_folds = 5,
                       alpha = 1/4,
                       gamma = 1, 
                       n_threads = 1)
rmss_coefs <- coef(cv.rmss_fit, 
                   h_ind = cv.rmss_fit$h_opt, t_ind = cv.rmss_fit$t_opt, u_ind = cv.rmss_fit$u_opt,
                   group_index = 1:cv.rmss_fit$n_models)
sens_rmss <- sum(which((rmss_coefs[-1]!=0)) <= p.active)/p.active
spec_rmss <- sum(which((rmss_coefs[-1]!=0)) <= p.active)/sum(rmss_coefs[-1]!=0)
rmss_preds <- predict(cv.rmss_fit, newx = x_test,
                      h_ind = cv.rmss_fit$h_opt, t_ind = cv.rmss_fit$t_opt, u_ind = cv.rmss_fit$u_opt,
                      group_index = 1:cv.rmss_fit$n_models,
                      dynamic = FALSE)
mean((y_test - rmss_preds)^2)/sigma^2
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
