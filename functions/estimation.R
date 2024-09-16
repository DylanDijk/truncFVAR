library(plyr)

#### Helper Functions ####
# Computes autocovariance without centering
Rcpp::cppFunction("
// Define sample covariance function without centering, with lag
NumericMatrix acf_no_center(NumericMatrix data, int lag = 0) {
  int n = data.nrow(); // Number of rows in the data matrix
  int p = data.ncol(); // Number of columns in the data matrix
  NumericMatrix cov_mat(p, p); // Covariance matrix
  
  // Loop through each pair of columns
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      double sum = 0.0; // Sum of product of elements
      // Compute sum of product of elements for each pair of rows with lag
      for (int k = 0; k < n - lag; k++) {
        sum += data(k+lag, i) * data(k, j);
      }
      // Compute sample covariance without centering with lag
      cov_mat(i, j) = sum / (n - lag);
    }
  }

  return cov_mat;
}")

# generates vector of median absolute deviation for each column of data matrix
mad_variables = function(data){
  mad_data = vector(length = ncol(data))
  for(i in 1:ncol(data)){
    mad_data[i] = mad(data[,i])
  }
  return(mad_data)
}

# given a matrix it creates a linear grid of taus from
# largest absolute value to median absolute value.
# In practice the function should be used on standardised data.
# Returns a list of two objects the grid of taus and then the quantiles
tau_grid_fun = function(data, n_steps){
  
  max_dat = max(abs(data))
  median_dat = median(abs(data))
  tau_grid = seq(median_dat, max_dat, length.out = n_steps)
  tau_grid_quant = sapply(tau_grid, function(x){ecdf(abs(data))(x)})
  
  return(list(tau_grid, tau_grid_quant))
}

tau_grid_stand_fun = function(data, n_steps, standardise){
  if(standardise){
    mad_data = mad_variables(data)
    tau_values = tau_grid_fun( t(t(data)/mad_data), n_steps = n_steps)
    mad_data_mat = matrix(mad_data, nrow = length(tau_values[[1]]), ncol = length(mad_data), byrow = T)
    tau_grid = mad_data_mat * tau_values[[1]]
  } else {
    # just a matrix of 1s no scaling
    mad_data_mat = matrix(1, nrow = n_steps, ncol = ncol(data), byrow = T)
    tau_values = tau_grid_fun(data, n_steps = n_steps)
    tau_grid = mad_data_mat * tau_values[[1]]
  }
  
  return(list(tau_values = tau_values, tau_grid = tau_grid))
}


# Computes autocovariance or covariance depending on lag argument, without centering.
# After truncating the data.
Rcpp::cppFunction("
NumericMatrix truncateAndComputeCovariance_lag(NumericMatrix data, NumericVector tau, int lag = 0) {
  int n = data.nrow(); // Number of rows in the data matrix
  int p = data.ncol(); // Number of columns in the data matrix
  NumericMatrix data_copy = Rcpp::clone(data); // Create a copy of the data matrix
  NumericMatrix covariance(p, p); // Matrix to store the covariance
  
  // Loop through each row of the data matrix
  for (int i = 0; i < n; i++) {
    // Loop through each column of the data matrix
    for (int j = 0; j < p; j++) {
      // Truncate the data element using the corresponding tau value
      if (std::abs(data_copy(i, j)) > tau[j]) {
        data_copy(i, j) = tau[j] * std::copysign(1.0, data_copy(i, j)); // Truncate to tau with sign
      }
    }
  }
  
  // Compute the covariance without estimating the sample mean
  if (lag > 0) {
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        double sum = 0.0;
        for (int k = lag; k < n; k++) {
          sum += data_copy(k, i) * data_copy(k-lag, j);
        }
        covariance(i, j) = sum / (n - lag);
      }
    }
  } else {
    for (int i = 0; i < p; i++) {
      for (int j = i; j < p; j++) { // Only need to compute the upper triangular part
        double sum = 0.0;
        for (int k = 0; k < n; k++) {
          sum += data_copy(k, i) * data_copy(k, j);
        }
        covariance(i, j) = sum / (n);
        covariance(j, i) = covariance(i, j); // Covariance matrix is symmetric
      }
    }
  }
  
  return covariance;
}")


#### Cross-Validation ####

# takes a single data matrix.
# creates grid of tau values on standardised data then rescales

# then data is split into two halves
 # for each tau value in the grid of taus 
 # the truncated cov is calculated on each half and so is the sample cov
 # then the frobenius norm of the difference of each estimate is calculated
 # then there are trim arguments to decided whether to trim large values in the difference matrix
# max = F, means that the MAX norm is NOT used in CV but the frobenius norm is used.
# If standardise is true, then taus are scaled by MAD of variables
cross_val = function(data, n_tau = 60, lag, trim = 1, trim_d = 1, trim_all = F, max, standardise = T){
  
  # get tau grid, standardised or not
  tau_l = tau_grid_stand_fun(data = data, n_steps = n_tau, standardise = standardise)
  # split data
  n = nrow(data)
  half_1 = data[1:(n/2),]
  half_2 = data[((n/2) + 1):n,]
  
  half_1_trunc = alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_1, lag = lag)
  half_2_trunc = alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_2, lag = lag)
  
  half_1_sample = acf_no_center(data = half_1, lag = lag)
  half_2_sample = acf_no_center(data = half_2, lag = lag)

  if(!max){ 
    trimmed_frob_norm =  function(x, mat_2){
      diff_mat = abs(x - mat_2)
      
      if(trim_all){
        if(trim < 1){
          diff_mat[diff_mat > quantile(diff_mat[upper.tri(diff_mat, diag = T)], trim)] = 0
        } 
        # Now looking to trim diagonal and off diagonal entries separately
      } else {
        
        if(trim < 1){
          if(lag == 0){
            diff_mat[(upper.tri(diff_mat) |  lower.tri(diff_mat))
                     & (diff_mat > quantile(diff_mat[upper.tri(diff_mat)], trim))] = 0
          } else {
            diff_mat[(upper.tri(diff_mat) |  lower.tri(diff_mat))
                     & (diff_mat > quantile(diff_mat[upper.tri(diff_mat) |  lower.tri(diff_mat)], trim))] = 0
          }
        }
      
        if(trim_d < 1){
          diag(diff_mat)[diag(diff_mat) > quantile(diag(diff_mat), trim_d)] = 0
        }
      }
      norm(diff_mat, "F")
    }
    
    half_1_trunc_err = lapply(half_1_trunc, trimmed_frob_norm, half_2_sample)
    half_2_trunc_err = lapply(half_2_trunc, trimmed_frob_norm, half_1_sample)
    
  } else {
    half_1_trunc_err = lapply(half_1_trunc, function(x){norm((x - half_2_sample), "M")})
    half_2_trunc_err = lapply(half_2_trunc, function(x){norm((x - half_1_sample), "M")})
  }
  
  
  scores_per_tau = 1/2 * (unlist(half_1_trunc_err) + unlist(half_2_trunc_err))
  names(scores_per_tau) = tau_l$tau_values[[1]]
  
  return(list(scores_per_tau, tau_l$tau_values[[2]]))
  
}


# calls the cross_val() function above and returns the tau vector that produce the smallest cross-validation error.
# then computes the error of the truncated estimate using the full dataset against the true covariance which we know
# with simulated data.

# is applied to a list of matrices
# cv_lag argument is to decide whether to use the autocovariances within the CV measure
# cov_error_comb argument is to decide whether to compute the error across lags as in the \mathcal{E}_{n, p} set.
cross_val_and_error = function(data, n_tau, lag, error = NULL, trim, trim_d, trim_all, max, standardise = T, true_cov, error_norm,
                               cv_lag=F, cov_error_comb = F){
  
  if(!cv_lag){
    scores_per_tau = lapply(data, cross_val, n_tau, lag, trim, trim_d, trim_all, max, standardise = standardise)
  }else{
    # default set to getting cv error using up to lag 1, if cv_lag = T
    scores_per_tau = lapply(data, cross_val_lag, n_tau, lag = 1, trim, trim_d, trim_all, max, standardise = standardise)
  }
  min_tau = lapply(scores_per_tau, function(x){as.numeric(names(which.min(x[[1]])))})
  min_tau_quant = lapply(scores_per_tau, function(x){x[[2]][which.min(x[[1]])]})
  
  if(is.null(error)){
    names(min_tau) = min_tau_quant
    return(unlist(min_tau))
  } else {
    if(standardise){
      mad_data = lapply(data, mad_variables)
      min_tau_vec = Map('*', mad_data, min_tau)
    } else {
      # not standardising so just have the same tau value for each variable
      min_tau_vec = lapply(min_tau, function(x) rep(x, ncol(data[[1]])))
    }
    
    if(cov_error_comb){
      # creating list of combined covariance and autocovariance matrices, by column binding
      truncated_cov = mapply(truncateAndComputeCovariance_lag, tau = min_tau_vec, data = data, lag = 0, SIMPLIFY = F)
      truncated_auto_cov = mapply(truncateAndComputeCovariance_lag, tau = min_tau_vec, data = data, lag = 1, SIMPLIFY = F)
      truncated_cov = mapply(function(m1, m2) cbind(m1, m2), truncated_cov, truncated_auto_cov, SIMPLIFY = FALSE)
      
    }else{
      truncated_cov = mapply(truncateAndComputeCovariance_lag, tau = min_tau_vec, data = data, lag = lag, SIMPLIFY = F)
    }
    
    if(is.function(error)){
      errors = lapply(truncated_cov, error)
    } else {
      errors = Map(function(x,y){norm(x-y, error_norm)}, truncated_cov, true_cov)
    }
    
    names(errors) = min_tau
  
    return(list((unlist(errors)), unlist(min_tau_quant)))
  }  
}


# takes a data matrix
# And returns a truncated data matrix where tau has been chosen by cross-validation.
# Previously, n_tau, lag, trim, trim_d, trim_all did not have default arguments set.
cross_val_and_trunc = function(data, n_tau = 60, lag = 0, trim = 1, trim_d = 1, trim_all = F, max = T, cv_lag= F, standardise = T){
  
  if(!cv_lag){
    tau_scores = cross_val(data, n_tau, lag, trim, trim_d, trim_all, max = max, standardise = standardise)
  }else{
    # default set to getting cv error using up to lag 1, if cv_lag = T
    tau_scores = cross_val_lag(data, n_tau, lag = 1, trim, trim_d, trim_all, max = max, standardise = standardise)
  }
  # tau_scores = cross_val(data, n_tau, lag, trim, trim_d, trim_all, max = max)
  min_tau = as.numeric(names(which.min(tau_scores[[1]])))
  
  if(standardise){
    mad_data = mad_variables(data)
    tau_s = mad_data * min_tau
  } else {
    tau_s = rep(1, ncol(data)) * min_tau
  }
  
  for(i in 1:nrow(data)){
    x = data[i,]
    data[i,] = ifelse(abs(x) > abs(tau_s), sign(x)*tau_s, x)
  }
  
  return(data)
  
}


##### Looking at max norm error including the acf for a lag, similar to error in equation 11 of fnets
# modified version of cross_val function
cross_val_lag = function(data, n_tau, lag, trim = 1, trim_d = 1, trim_all = F, max = T, standardise = standardise){
  
  # get tau grid
  tau_l = tau_grid_stand_fun(data = data, n_steps = n_tau, standardise = standardise)
  
  # split data
  n = nrow(data)
  half_1 = data[1:(n/2),]
  half_2 = data[((n/2) + 1):n,]
  
  scores_per_tau = vector(mode = "list", length = lag+1)
  
  for(i in 0:lag){
    half_1_trunc = alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_1, lag = i)
    half_2_trunc = alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_2, lag = i)
  
    half_1_sample = acf_no_center(data = half_1, lag = i)
    half_2_sample = acf_no_center(data = half_2, lag = i)
    
    if(max == T){
      half_1_trunc_err = lapply(half_1_trunc, function(x){norm((x - half_2_sample), "M")})
      half_2_trunc_err = lapply(half_2_trunc, function(x){norm((x - half_1_sample), "M")})
    }else{
      half_1_trunc_err = lapply(half_1_trunc, function(x){norm((x - half_2_sample), "F")})
      half_2_trunc_err = lapply(half_2_trunc, function(x){norm((x - half_1_sample), "F")})
    }
      
    scores_per_tau[[i+1]] = 1/2 * (unlist(half_1_trunc_err) + unlist(half_2_trunc_err))
  }
  
  scores_per_tau = do.call(pmax, scores_per_tau)
  
  names(scores_per_tau) = tau_l$tau_values[[1]]
  
  return(list(scores_per_tau, tau_l$tau_values[[2]]))
  
}




#### Computing the covariance of a VAR 1 model ####
# Lutkephol book equation 2.1.32 gives formula for calculating cov matrix.
# This requires estimating the inverse, which becomes infeasible for large dimension.
# A paper by Barone 1987 https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1467-9892.1987.tb00426.x
# mentions an iterative algorithm from a book of Anderson and Moore 1979.

# I now provide a function that calculates this 
library(Matrix)
cov_of_var = function(tol = 2.220446e-16, A, print = F){
  p = ncol(A)
  M_old = A
  N_old = Diagonal(p)
  converged = F
  
  while(!converged) {
    N_new = M_old %*% N_old %*% t(M_old) + N_old
    M_new = M_old %*% M_old
    change = norm(N_new - N_old, "M")
    if(print == T){
      print(change)
    }
    M_old = M_new
    N_old = N_new
    converged = ifelse(change < tol, T, F)
  }
  return(N_old)
}

# Now computing the covariance by the exact equation which is feasible for coefficient matrices of small dimension.
cov_of_var_manual = function(A){
  
  A = Matrix(A, sparse = T)
  p = ncol(A)
  inv = Matrix::solve((Diagonal(p^2) - (kronecker(A,A))))
  inv = Matrix(inv, sparse = T)
  cov_var = matrix(round(inv %*% sparseVector(rep(1,p), seq(0,(p^2-1),p) + 1:p, p^2),3), p, p)
  return(cov_var)
  
}



######## static.pca ######## 
# fnets version uses data in p x n format as input
m.static.pca = function(xx, q, mm = 1){
  n <- dim(xx)[1]
  p <- dim(xx)[2]
  covx <- t(xx) %*% xx / n
  eig <- eigen(covx, symmetric = TRUE)
  
  if(q >= 1){
    proj <- eig$vectors[, 1:q, drop = FALSE] %*% t(eig$vectors[, 1:q, drop = FALSE])
    lam <- eig$vectors[, 1:q, drop = FALSE] %*% diag(eig$values[1:q])
    f <- xx %*% (eig$vectors[, 1:q, drop = FALSE]) %*% diag(1/(eig$values[1:q]))
  } else lam <- f <- NULL
  
  Gamma_c <- Gamma_x <- array(0, dim = c(p, p, 2 * mm + 1))
  for(h in 0:mm) {
    ifelse(h == 0, Gamma_x[,, h + 1] <- covx, Gamma_x[,, h + 1] <- t(xx[1:(n - h), ]) %*% xx[(1:(n - h)) + h, ] / n)
    if(q >= 1) Gamma_c[,, h + 1] <- proj %*% Gamma_x[,, h + 1] %*% proj
    if(h != 0) {
      Gamma_x[,, 2 * mm + 1 - h + 1] <- t(Gamma_x[,, h + 1])
      Gamma_c[,, 2 * mm + 1 - h + 1] <- t(Gamma_c[,, h + 1])
    }
  }
  acv <- list(Gamma_x = Gamma_x, Gamma_c = Gamma_c, Gamma_i = Gamma_x - Gamma_c)
  
  out <- list(q = q, lam = lam, f = f, acv = acv)
  return(out)
  
}


glm_pca_fun = function(x, q){
  if(q > 0){
    mpca = m.static.pca(x, q = q)
    xx_c = mpca$f %*% t(mpca$lam)
    xx_i = x - xx_c
    glmn = sparsevar::fitVAR(xx_i, p = 1, method = "cv", parallel = TRUE, ncores = 4)
  } else {
    glmn = sparsevar::fitVAR(x, p = 1)
  }  
  glmn_A = glmn$A[[1]]
  return(list(A = glmn_A, glmn = glmn, common = xx_c))
}



