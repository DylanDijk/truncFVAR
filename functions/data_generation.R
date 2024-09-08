# Generate data from VAR(1) process, with innovations from standardised t distributions
# Either p is a sequence or n.
# A list of length of the sequence is generated, where each element contains nsim number of simulations

#individual functions for generating errors from distributions
stand_t = function(n, df){
  sd_df = sqrt(df/(df-2))
  return(rt(n = n, df = df)/sd_df)
}

std_log_norm = function(meanlog = 0, sdlog = 1, n =1){
  return((rlnorm(n=n, meanlog = meanlog, sdlog = sdlog) - exp(1/2))/sqrt(exp(2) - exp(1)))
}

add_innov = function(data, innov_dist, covstr, innov_df, p){

  if(innov_dist == "t"){
    innov = stand_t(n = p, df = innov_df)
  }
  if(innov_dist == "log"){
    innov = (std_log_norm(n = p))
  }
  if(innov_dist == "gauss"){
    innov = rnorm(n = p)
  }
  
  if(!is.null(covstr)){
    innov = covstr %*% innov
  }
  
  data = data + innov
  return(data)
}

#####
# single function that takes innov_dist as argument
VAR_1_data_ind = function(nsim, n_p, A_coeff, innov_df, covstr = NULL, innov_dist = "t"){
  
  sim_VAR_data = vector(mode = "list", length = nrow(n_p))
  
  list_names <- character(length = nrow(n_p))
  
  if(!is.list(A_coeff)){
    warning("A_coeff is not a list. Using a random generator function. Default is rand_D_mat.")
    # for example when using D matrix as described in report for VAR process to generate factor process.
    # The function to generate the coefficient is called later on in the function at A = A_coeff()
  }
  
  # Loop through each row of the matrix and create the strings
  for (i in 1:nrow(n_p)) {
    list_names[i] <- paste("(", n_p[i, "n"], ",", n_p[i, "p"], ")", sep = "")
  }
  
  names(sim_VAR_data) = list_names
  
  for(i in 1:nrow(n_p)){
    # sim_VAR_data[[list_names[i]]] = vector(mode = "list", length = nsim)
    sim_VAR_data[[i]] = vector(mode = "list", length = nsim)
    attr(sim_VAR_data[[i]], "A") = vector(mode = "list", length = nsim)
    
    n = n_p[i, "n"]
    p = n_p[i, "p"]
    
    for(k in 1:nsim){
      
      if(is.list(A_coeff)){
        A = A_coeff[[i]]
      }else{
        # case where we are randomly sampling the coeff matrix, e.g for the factor model simulations
        A = A_coeff(p = p)
        attr(sim_VAR_data[[i]], "A")[[k]] = A
      }
      
      data = matrix(nrow = n, ncol = p)
      data[1,] = add_innov(data = rep(0, ncol(data)), innov_dist = innov_dist, covstr = covstr[[i]], p = p, innov_df = innov_df)
      for(j in 2:n){
        
        data[j,] = A %*% data[j-1, ]
        data[j,] = add_innov(data = data[j,], innov_dist = innov_dist, covstr = covstr[[i]], p = p, innov_df = innov_df)
        
      }
      # sim_VAR_data[[list_names[i]]][[k]] = data
      sim_VAR_data[[i]][[k]] = data
    }
  }
  return(sim_VAR_data)
}



# Generate standardised independent data from t distribution.
# Each vector independent, and each coordinate of vector independent.
# n_p argument of function is a matrix determining the values of n and p that should be used in the simulation
# For example, n_p_values = cbind(n = c(200), p = seq(50, 300, 50))
ind_data_ind_t = function(nsim, n_p, df, covstr = NULL){
  
  sd_df = sqrt(df/(df-2))
  
  sim_t_data = vector(mode = "list", length = nrow(n_p))
  
  list_names <- character(length = nrow(n_p))
  
  # Loop through each row of the matrix and create the strings
  for (i in 1:nrow(n_p)) {
    list_names[i] <- paste("(", n_p[i, "n"], ",", n_p[i, "p"], ")", sep = "")
  }
  
  names(sim_t_data) = list_names
  
  if(is.null(covstr)){
    for(i in 1:nrow(n_p)){
      sim_t_data[[list_names[i]]] = vector(mode = "list", length = nsim)
      
      n = n_p[i, "n"]
      p = n_p[i, "p"]
      
      for(j in 1:nsim){
        data = rt(n = n*p, df = df)/(sd_df)
        sim_t_data[[list_names[i]]][[j]] = matrix(data, nrow = n, ncol = p)
      }
    }
  } else {
    for(i in 1:nrow(n_p)){
      sim_t_data[[list_names[i]]] = vector(mode = "list", length = nsim)
      
      n = n_p[i, "n"]
      p = n_p[i, "p"]
      
      for(j in 1:nsim){
        data = rt(n = n*p, df = df)/(sd_df)
        sim_t_data[[list_names[i]]][[j]] = matrix(data, nrow = n, ncol = p) %*% covstr[[i]]
      }
    }
  }
  return(sim_t_data)
}





# Generate independent data from standard Gaussian distribution.
# Each vector independent, and each coordinate of vector independent.
# covstr argument is so that a covariance structure can be added to the data.
# this takes the form of a list covstr = vector(mode = "list", length = nrow(n_p_values))
ind_data_ind_gauss = function(nsim, n_p, covstr = NULL){
  sim_gauss_data = vector(mode = "list", length = nrow(n_p))
  
  list_names <- character(length = nrow(n_p))
  
  # Loop through each row of the matrix and create the strings
  for (i in 1:nrow(n_p)) {
    list_names[i] <- paste("(", n_p[i, "n"], ",", n_p[i, "p"], ")", sep = "")
  }
  
  names(sim_gauss_data) = list_names
  
  if(is.null(covstr)){
    for(i in 1:nrow(n_p)){
      sim_gauss_data[[list_names[i]]] = vector(mode = "list", length = nsim)
      
      n = n_p[i, "n"]
      p = n_p[i, "p"]
      
      for(j in 1:nsim){
        data = rnorm(n = n*p)
        sim_gauss_data[[list_names[i]]][[j]] = matrix(data, nrow = n, ncol = p)
      }
    }
  } else {
    for(i in 1:nrow(n_p)){
      sim_gauss_data[[list_names[i]]] = vector(mode = "list", length = nsim)
      
      n = n_p[i, "n"]
      p = n_p[i, "p"]
      
      for(j in 1:nsim){
        data = rnorm(n = n*p)
        sim_gauss_data[[list_names[i]]][[j]] = matrix(data, nrow = n, ncol = p) %*% covstr[[i]]
      }
    }
  }
    
  return(sim_gauss_data)
}

ind_data_ind_lognorm = function(nsim, n_p, covstr = NULL){
  sim_gauss_data = vector(mode = "list", length = nrow(n_p))
  
  list_names <- character(length = nrow(n_p))
  
  # Loop through each row of the matrix and create the strings
  for (i in 1:nrow(n_p)) {
    list_names[i] <- paste("(", n_p[i, "n"], ",", n_p[i, "p"], ")", sep = "")
  }
  
  names(sim_gauss_data) = list_names
  
  if(is.null(covstr)){
    for(i in 1:nrow(n_p)){
      sim_gauss_data[[list_names[i]]] = vector(mode = "list", length = nsim)
      
      n = n_p[i, "n"]
      p = n_p[i, "p"]
      
      for(j in 1:nsim){
        data = std_log_norm(n = n*p)
        sim_gauss_data[[list_names[i]]][[j]] = matrix(data, nrow = n, ncol = p)
      }
    }
  } else {
    for(i in 1:nrow(n_p)){
      sim_gauss_data[[list_names[i]]] = vector(mode = "list", length = nsim)
      
      n = n_p[i, "n"]
      p = n_p[i, "p"]
      
      for(j in 1:nsim){
        data = std_log_norm(n = n*p)
        sim_gauss_data[[list_names[i]]][[j]] = matrix(data, nrow = n, ncol = p) %*% covstr[[i]]
      }
    }
  }
  
  return(sim_gauss_data)
}




# Generates list of banded covariance matrices based on list input n_p
A_coeff_banded = function(n_p){
  A_coeff = vector(mode = "list", length = nrow(n_p))
  for(k in 1:nrow(n_p)){
    
    names(A_coeff)[k] <- paste("(", n_p[k, "n"], ",", n_p[k, "p"], ")", sep = "")
    
    p = n_p[k, "p"]
    
    mat = matrix(0, nrow = p, ncol = p) # initialize matrix with all elements set to 0
    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) {
          mat[i,j] <- 0.5 # set diagonal elements to 0.5
        } else if (i-j == 1) {
          mat[i,j] <- 0.4 # set elements above the diagonal to 0.4
        } else if (i-j == -1) {
          mat[i,j] <- -0.4 # set elements below the diagonal to -0.4
        }
      }
      A_coeff[[k]] = mat
    }
  }
  return(A_coeff)
}

# generates a (random) renyi graph coefficient matrix
renyi_coeff = function(p){
  n = 2
  A = fnets::sim.var(n,p)$A
  return(A)
}



#### Functions used for factor plus var process simulations   ####

# D, coeff matrix for VAR model generating factors
# Barigozzi, Cho, and Owens (2022, Section D, C2)
rand_D_mat = function(p){
  r = p
  D_coeff = matrix(runif(r*r,0,0.3), r, r)
  diag(D_coeff) = runif(r, 0.5, 0.8)
  D_coeff = 0.7 * D_coeff
  # check largest eigenvalue of D_coeff
  eigen_res = eigen(D_coeff)
  max(Mod(eigen_res$values))
  if(max(Mod(eigen_res$values)) > 1){
    D_coeff = D_coeff / max(Mod(eigen_res$values))
  }
  return(D_coeff)
}




# this function generates data from the factor model with VAR errors
# requires A, Lambda, and D.
# If norm_loading is T then loading entries are normalised
factor_mod_data_VAR_idio = function(dist = "t", innov_df = 2.1, indo_A = A_coeff, n_p, fac_var_ret = F, norm_loading = F){
  
  var_factor = VAR_1_data_ind(nsim = nsim, n_p = n_r, A_coeff = rand_D_mat, innov_df = innov_df, innov_dist = dist)
  var_indo = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = indo_A, innov_df = innov_df, innov_dist = dist)
  
  n_p_names = names(A_coeff_banded(n_p))
  names(var_factor) = n_p_names
  
  if(norm_loading){
    A_norm = vector(mode = "list", length = nrow(n_p)); names(A_norm) = names(A_coeff_banded(n_p))
    
    for(i in 1:nrow(n_p)){
      A_norm[[i]] = matrix(nrow = nsim, ncol = n_p[i,"p"])
    }
    
    # sqrt(
    # diag(cov_of_var(A = attributes(var_indo$`(200,50)`)$A[[104]])) /
    # diag(t(lambda_list$`(200,50)`[[104]]) %*% cov_of_var(A = attributes(var_factor[[3]])$A[[104]]) %*% lambda_list$`(200,50)`[[104]])
    # )
    
    A_norm_fun = function(indo, fac, A_norm, lambda_list, A_coeff){
      for(i in 1:nsim){
        if(is.function(indo_A)){
          A_norm[i,] =
          sqrt(
            diag(cov_of_var(A = attributes(indo)$A[[i]])) /
            diag(t(lambda_list[[i]]) %*% cov_of_var(A = attributes(fac)$A[[i]]) %*% lambda_list[[i]])
          )
        }else{
          A_norm[i,] =
            sqrt(
              diag(cov_of_var(A = A_coeff)) /
                diag(t(lambda_list[[i]]) %*% cov_of_var(A = attributes(fac)$A[[i]]) %*% lambda_list[[i]])
            )
        }  
      }
      return(A_norm)
    }
    # browser()
    for(i in 1:nrow(n_p)){
      A_norm[[i]] = A_norm_fun(indo = var_indo[[i]], fac = var_factor[[i]], A_norm = A_norm[[i]], lambda_list = lambda_list[[i]], A_coeff = indo_A[[i]])
    }
    
  }
  
  
  data_comb = vector(mode = "list", length = length(n_p_names))
  names(data_comb) = n_p_names
  
  for(i in 1:nrow(n_p)){
    data_comb[[n_p_names[i]]] = vector(mode = "list", length = nsim)
    
    for(k in 1:nsim){
      if(norm_loading){
        data = sweep(var_factor[[i]][[k]] %*% lambda_list[[i]][[k]], 2, A_norm[[i]][k,], `*`) + var_indo[[i]][[k]]
      } else {  
        data = (var_factor[[i]][[k]] %*% lambda_list[[i]][[k]]) + var_indo[[i]][[k]]
      }
      data_comb[[n_p_names[i]]][[k]] = data
    }
  }
  
  attr(data_comb, "A_fac") = lapply(var_factor, function(x){attr(x, "A")})
  if(fac_var_ret){
    attr(data_comb, "fac_var") = var_factor
  }
  if(is.function(indo_A)){
    attr(data_comb, "A_var") = lapply(var_indo, function(x){attr(x, "A")})
  }
  if(norm_loading){
    attr(data_comb, "A_norm") = A_norm
  }
  
  return(data_comb)
  
}

##### For cases when I just want to quickly generate a single simulations of a factor plus var dataset.
# I have the single function below.
# Simulated as described in (F1) in report
  # VAR coefficient matrix (A in report) is banded 
  # lambda all standard normal
  # D matrix, as described
fac_var_dat = function(n,p,r, dist = "gauss", innov_df, A = "banded"){
  n_p = cbind(n = n, p = p)
  nsim = 1
  
  if(r > 0){
    n_r = cbind(n = n, p = r)
    lambda_list = vector(mode = "list", length = 1)
    names(lambda_list) <-p
    lambda_list[[1]][[1]] = matrix(rnorm(r*p), nrow = r,  ncol = p)
  }
  if(A == "banded"){
    A_coeff = A_coeff_banded(n_p)
  } 
  if(A == "renyi"){
    A_coeff_name = names(A_coeff_banded(n_p))
    A_coeff = list()
    A_coeff[[A_coeff_name]] = fnets::sim.var(n,p)$A
  }
  var_indo = VAR_1_data_ind(nsim = nsim, n_p = n_p, A_coeff = A_coeff, innov_df = innov_df, innov_dist = dist)
  
  if(r > 0){
    var_factor = VAR_1_data_ind(nsim = nsim, n_p = n_r, A_coeff = rand_D_mat, innov_df = innov_df, innov_dist = dist)
    data = var_factor[[1]][[1]] %*% lambda_list[[1]][[1]] + var_indo[[1]][[1]]
    attr(data, "A") = A_coeff
    attr(data, "r") = r
    return(data)
  } else {
    data = var_indo[[1]][[1]]
    attr(data, "A") = A_coeff
    attr(data, "r") = r
    return(data)
  }  
}










