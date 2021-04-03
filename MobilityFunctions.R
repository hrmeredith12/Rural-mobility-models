## Functions called in "Modeling Rural Mobility in Burkina Faso" notebook

## Create matrix with origin as rows and destination as columns for various variables (trips, trip types)
make.OD.matrix <- function(datafile, origin, destination, var.of.interest, IDs){
  
  monthly.trips <- datafile%>%
    pivot_wider(id_cols = c(origin, destination),
                names_from = destination,
                values_from = var.of.interest)
  monthly.trips <- as.matrix(monthly.trips[ ,-1])
  row.names(monthly.trips) <- IDs
  
  return(monthly.trips)
  
}

##### Gravity models #####

#### Basic Gravity model with power distance function ####
fit.basic.gravity.model.pwr <- function(M, D, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function(M, D, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                           n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                           parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <-"   
      model {
        ## Poisson likelihood
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            M[i,j] ~ dpois(c[i,j])  
          }
        }
        
        ## Gravity model
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            c[i,j] <- ifelse(i == j, 1e-06,
                             exp(log(theta) + (omega_1*log(N_dest[i]) + omega_2*log(N_orig[j]) - log( f_d[i,j] )))
            )
            f_d[i,j] <- D[i,j]^gamma
          }
        }
        
        ## Priors
        theta ~ dgamma(prior_theta[1], prior_theta[2])
        omega_1 ~ dgamma(prior_omega_1[1], prior_omega_1[2])
        omega_2 ~ dgamma(prior_omega_2[1], prior_omega_2[2])
        gamma ~ dgamma(prior_gamma[1], prior_gamma[2])
      }"
    
    params <- c("omega_1", "omega_2", "theta","gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
    
  }   
    
    output <- fit.gravity.model(M = monthly.trips,
                                D = distance,
                                N = N,
                                n_chain = n_chain,
                                n_burn = n_burn,
                                n_samp = n_samp,
                                n_thin = n_thin,
                                prior = NULL,
                                DIC = TRUE,
                                parallel = TRUE)
    
    parameters <- mobility::summarize_mobility(output)
    return(parameters)
}

#### Basic Gravity model with exponential distance function ####
fit.basic.gravity.model.exp <- function(M, D, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function(M, D, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                                n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                                parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <-"   
      model {
        ## Poisson likelihood
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            M[i,j] ~ dpois(c[i,j])  
          }
        }
        
        ## Gravity model
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            c[i,j] <- ifelse(i == j, 1e-06,
                             exp(log(theta) + (omega_1*log(N_dest[i]) + omega_2*log(N_orig[j]) - log( f_d[i,j] )))
            )
            f_d[i,j] <- exp(D[i,j]/gamma)
          }
        }
        
        ## Priors
        theta ~ dgamma(prior_theta[1], prior_theta[2])
        omega_1 ~ dgamma(prior_omega_1[1], prior_omega_1[2])
        omega_2 ~ dgamma(prior_omega_2[1], prior_omega_2[2])
        gamma ~ dgamma(prior_gamma[1], prior_gamma[2])
      }"
    
    params <- c("omega_1", "omega_2", "theta","gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
    
  }   
  
  output <- fit.gravity.model(M = monthly.trips,
                              D = distance,
                              N = N,
                              n_chain = n_chain,
                              n_burn = n_burn,
                              n_samp = n_samp,
                              n_thin = n_thin,
                              prior = NULL,
                              DIC = TRUE,
                              parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}


#### Regional Gravity model with power distance function ####
fit.regional.gravity.model.pwr <- function(M, D, TT, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function (M, D, TT, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                                 n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                                 parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      TT = TT,
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <- "  
  model {
      # Poisson likelihood
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
  
           M[i,j] ~ dpois(c[i,j])  
         }
      }
  
      # Gravity model
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
          c[i,j] <- ifelse(
            i == j, 1e-06,
             ifelse(TT[i, j] == 1, 
                   exp(log(theta) + (omega_1[1] * log(N_dest[i]) + omega_2[1] * log(N_orig[j]) - gamma[1] * log(D[i,j]))),
                   exp(log(theta) + (omega_1[2] * log(N_dest[i]) + omega_2[2] * log(N_orig[j]) - gamma[2] * log(D[i,j])))
            )
         )
        }
      }
  
      # Priors
      theta ~ dgamma(prior_theta[1], prior_theta[2])
  
      for (k in 1:2) { #[1] = Intra-regional; [2] = Inter-regional
        omega_1[k] ~ dgamma(prior_omega_1[1], prior_omega_1[2])
        omega_2[k] ~ dgamma(prior_omega_2[1], prior_omega_2[2])
        gamma[k] ~ dgamma(prior_gamma[1], prior_gamma[2])
      }
    }"
    params <- c("omega_1", "omega_2", "theta", "gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
  }
  
  output <- fit.gravity.model(M = monthly.trips,
                              D = distance,
                              TT = TT, 
                              N = N,
                              n_chain = n_chain,
                              n_burn = n_burn,
                              n_samp = n_samp,
                              n_thin = n_thin,
                              prior = NULL,
                              DIC = TRUE,
                              parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}

#### Regional Gravity model with exponential distance function ####
fit.regional.gravity.model.exp <- function(M, D, TT, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function (M, D, TT, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                                 n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                                 parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      TT = TT,
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <-"   
      model {
        # Poisson likelihood
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            
            M[i,j] ~ dpois(c[i,j])  
          }
        }
        
        # Gravity model
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            c[i,j] <- ifelse(
              i == j,1e-06,
              ifelse(TT[i, j] == 1, 
                     exp(log(theta) + (omega_1[1] * log(N_dest[i]) + omega_2[1] * log(N_orig[j]) - (D[i,j])/gamma[1])),
                     exp(log(theta) + (omega_1[2] * log(N_dest[i]) + omega_2[2] * log(N_orig[j]) - (D[i,j])/gamma[2]))))
          }
        }
        
        # Priors
        theta ~ dgamma(prior_theta[1], prior_theta[2])
        
        for (k in 1:2) { #[1] = Intra-regional; [2] = Inter-regional
          omega_1[k] ~ dgamma(prior_omega_1[1], prior_omega_1[2])
          omega_2[k] ~ dgamma(prior_omega_2[1], prior_omega_2[2])
          gamma[k] ~ dgamma(prior_gamma[1], prior_gamma[2])
        }
      }"
    params <- c("omega_1", "omega_2", "theta", "gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
  }
  
  output <- fit.gravity.model(M = monthly.trips,
                              D = distance,
                              TT = TT, 
                              N = N,
                              n_chain = n_chain,
                              n_burn = n_burn,
                              n_samp = n_samp,
                              n_thin = n_thin,
                              prior = NULL,
                              DIC = TRUE,
                              parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}

#### Urbanicity Gravity model with power distance function ####
fit.urban.gravity.model.pwr <- function(M, D, TT, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function (M, D, TT, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                                 n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                                 parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      TT = TT,
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <-"   
      model {
        # Poisson likelihood
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            M[i,j] ~ dpois(c[i,j])  
          }
        }
        
        # Gravity model
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            c[i,j] <- ifelse(i == j, 1e-06,            
                             ifelse(TT[i,j] == 1, exp(log(theta) + (omega_1[1] * log(N_dest[i]) + omega_2[1] * log(N_orig[j]) - gamma[1] * log(D[i,j]))),
                             ifelse(TT[i,j] == 2, exp(log(theta) + (omega_1[2] * log(N_dest[i]) + omega_2[2] * log(N_orig[j]) - gamma[2] * log(D[i,j]))),
                             ifelse(TT[i,j] == 3, exp(log(theta) + (omega_1[3] * log(N_dest[i]) + omega_2[3] * log(N_orig[j]) - gamma[3] * log(D[i,j]))),
                                                  exp(log(theta) + (omega_1[4] * log(N_dest[i]) + omega_2[4] * log(N_orig[j]) - gamma[4] * log(D[i,j])))))))
          }
        }
        
        # Priors
        theta ~ dgamma(prior_theta[1], prior_theta[2])
        
        for (k in 1:4) {
          omega_1[k] ~ dgamma(prior_omega_1[1], prior_omega_1[2])
          omega_2[k] ~ dgamma(prior_omega_2[1], prior_omega_2[2])
          gamma[k] ~ dgamma(prior_gamma[1], prior_gamma[2])
        }
      }"
    
    params <- c("omega_1", "omega_2", "theta", "gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
  }
  
  output <- fit.gravity.model(M = monthly.trips,
                              D = distance,
                              TT = TT, 
                              N = N,
                              n_chain = n_chain,
                              n_burn = n_burn,
                              n_samp = n_samp,
                              n_thin = n_thin,
                              prior = NULL,
                              DIC = TRUE,
                              parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}

#### Urbanicity Gravity model with exponential distance function ####
fit.urban.gravity.model.exp <- function(M, D, TT, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function (M, D, TT, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                           n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                           parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      TT = TT,
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <-"   
      model {
        # Poisson likelihood
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            M[i,j] ~ dpois(c[i,j])  
          }
        }
        
        # Gravity model
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
            c[i,j] <- ifelse(i == j, 1e-06,            
                             ifelse(TT[i,j] == 1, exp(log(theta) + (omega_1[1] * log(N_dest[i]) + omega_2[1] * log(N_orig[j]) - (D[i,j])/gamma[1])),
                             ifelse(TT[i,j] == 2, exp(log(theta) + (omega_1[2] * log(N_dest[i]) + omega_2[2] * log(N_orig[j]) - (D[i,j])/gamma[2])),
                             ifelse(TT[i,j] == 3, exp(log(theta) + (omega_1[3] * log(N_dest[i]) + omega_2[3] * log(N_orig[j]) - (D[i,j])/gamma[3])),
                                                  exp(log(theta) + (omega_1[4] * log(N_dest[i]) + omega_2[4] * log(N_orig[j]) - (D[i,j])/gamma[4]))))))
          }
        }
        
        # Priors
        theta ~ dgamma(prior_theta[1], prior_theta[2])
        
        for (k in 1:4) {
          omega_1[k] ~ dgamma(prior_omega_1[1], prior_omega_1[2])
          omega_2[k] ~ dgamma(prior_omega_2[1], prior_omega_2[2])
          gamma[k] ~ dgamma(prior_gamma[1], prior_gamma[2])
        }
      }"
    
    params <- c("omega_1", "omega_2", "theta","gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
  }
  
  output <- fit.gravity.model(M = monthly.trips,
                              D = distance,
                              TT = TT, 
                              N = N,
                              n_chain = n_chain,
                              n_burn = n_burn,
                              n_samp = n_samp,
                              n_thin = n_thin,
                              prior = NULL,
                              DIC = TRUE,
                              parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}

#### Regional-Urbanicity model with power distance function ####
fit.reg.urb.gravity.model.pwr <- function(M, D, TT, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function (M, D, TT, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                                 n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                                 parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      TT = TT,
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <- "  
  model {
  
      # Poisson likelihood
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
           M[i,j] ~ dpois(c[i,j])  
          }
      }
  
      # Gravity model
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
            c[i,j] <- ifelse(TT[i, j] == 0, 1e-06, 
                    ifelse(TT[i, j] == 1, exp(log(theta) + (omega_1[1] * log(N_dest[i]) + omega_2[1] * log(N_orig[j]) - gamma[1] * log(D[i,j]))),
                    ifelse(TT[i, j] == 2, exp(log(theta) + (omega_1[2] * log(N_dest[i]) + omega_2[2] * log(N_orig[j]) - gamma[2] * log(D[i,j]))),
                    ifelse(TT[i, j] == 3, exp(log(theta) + (omega_1[3] * log(N_dest[i]) + omega_2[3] * log(N_orig[j]) - gamma[3] * log(D[i,j]))),
                    ifelse(TT[i, j] == 4, exp(log(theta) + (omega_1[4] * log(N_dest[i]) + omega_2[4] * log(N_orig[j]) - gamma[4] * log(D[i,j]))),
                    ifelse(TT[i, j] == 5, exp(log(theta) + (omega_1[5] * log(N_dest[i]) + omega_2[5] * log(N_orig[j]) - gamma[5] * log(D[i,j]))),
                    ifelse(TT[i, j] == 6, exp(log(theta) + (omega_1[6] * log(N_dest[i]) + omega_2[6] * log(N_orig[j]) - gamma[6] * log(D[i,j]))),
                    ifelse(TT[i, j] == 7, exp(log(theta) + (omega_1[7] * log(N_dest[i]) + omega_2[7] * log(N_orig[j]) - gamma[7] * log(D[i,j]))),
                                         exp(log(theta) + (omega_1[8] * log(N_dest[i]) + omega_2[8] * log(N_orig[j]) - gamma[8] * log(D[i,j])))
                    ))))))))
          }
      }
  
      # Priors
      theta ~ dgamma(prior_theta[1], prior_theta[2])
  
      for (k in 1:8) { # note [0] = 'stay', [1]'IN-rural-rural', [2] = 'OUT-rural-rural', [3] = 'IN-rural-urban', [4] = 'OUT-rural-urban', [5] = 'IN-urban-rural', [6] = 'OUT-urban-rural', [7] = 'IN-urban-urban', [8] = 'OUT-urban-urban' 
        omega_1[k] ~ dgamma(prior_omega_1[1], prior_omega_1[2])
        omega_2[k] ~ dgamma(prior_omega_2[1], prior_omega_2[2])
        gamma[k] ~ dgamma(prior_gamma[1], prior_gamma[2])
      }
  
    }"
    params <- c("omega_1", "omega_2", "theta", "gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
  }
  output <- fit.gravity.model(M = monthly.trips,
                              D = distance,
                              TT = TT, 
                              N = N,
                              n_chain = n_chain,
                              n_burn = n_burn,
                              n_samp = n_samp,
                              n_thin = n_thin,
                              prior = NULL,
                              DIC = TRUE,
                              parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}

#### Regional-Urbanicity model with exponential distance function ####

fit.reg.urb.gravity.model.exp <- function(M, D, TT, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.gravity.model <- function (M, D, TT, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                                 n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                                 parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    message(paste("::Fitting gravity model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior, omega_1 = null_prior,
                    omega_2 = null_prior, gamma = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      D = D, 
                      TT = TT,
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta, 
                      prior_omega_1 = prior$omega_1,
                      prior_omega_2 = prior$omega_2, 
                      prior_gamma = prior$gamma)
    
    jags_model <- "  
     model {
    
        # Poisson likelihood
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
    
             M[i,j] ~ dpois(c[i,j])  
          }
    
        }
    
        # Gravity model
        for (i in 1:length(N_orig)) {
          for (j in 1:length(N_dest)) {
    
            
            c[i,j] <- ifelse(TT[i, j] == 0, 1e-06, 
                      ifelse(TT[i, j] == 1, exp(log(theta) + (omega_1[1] * log(N_dest[i]) + omega_2[1] * log(N_orig[j]) - (D[i,j])/gamma[1])),
                      ifelse(TT[i, j] == 2, exp(log(theta) + (omega_1[2] * log(N_dest[i]) + omega_2[2] * log(N_orig[j]) - (D[i,j])/gamma[2])),
                      ifelse(TT[i, j] == 3, exp(log(theta) + (omega_1[3] * log(N_dest[i]) + omega_2[3] * log(N_orig[j]) - (D[i,j])/gamma[3])),
                      ifelse(TT[i, j] == 4, exp(log(theta) + (omega_1[4] * log(N_dest[i]) + omega_2[4] * log(N_orig[j]) - (D[i,j])/gamma[4])),
                      ifelse(TT[i, j] == 5, exp(log(theta) + (omega_1[5] * log(N_dest[i]) + omega_2[5] * log(N_orig[j]) - (D[i,j])/gamma[5])),
                      ifelse(TT[i, j] == 6, exp(log(theta) + (omega_1[6] * log(N_dest[i]) + omega_2[6] * log(N_orig[j]) - (D[i,j])/gamma[6])),
                      ifelse(TT[i, j] == 7, exp(log(theta) + (omega_1[7] * log(N_dest[i]) + omega_2[7] * log(N_orig[j]) - (D[i,j])/gamma[7])),
                                           exp(log(theta) + (omega_1[8] * log(N_dest[i]) + omega_2[8] * log(N_orig[j]) - (D[i,j])/gamma[8]))
                      ))))))))
      
          }
        }
    
        # Priors
        theta ~ dgamma(prior_theta[1], prior_theta[2])
    
        for (k in 1:8) { # note [0] = 'stay',[1] = 'IN-rural-rural', [2] = 'OUT-rural-rural', [3] = 'IN-rural-urban', [4] = 'OUT-rural-urban', [5] = 'IN-urban-rural', [6] = 'OUT-urban-rural', [7] = 'IN-urban-urban', [8] = 'OUT-urban-urban' 
          omega_1[k] ~ dgamma(prior_omega_1[1], prior_omega_1[2])
          omega_2[k] ~ dgamma(prior_omega_2[1], prior_omega_2[2])
          gamma[k] ~ dgamma(prior_gamma[1], prior_gamma[2])
        }
    
      }"
    params <- c("omega_1", "omega_2", "theta", "gamma")
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
  }
  output <- fit.gravity.model(M = monthly.trips,
                              D = distance,
                              TT = TT, 
                              N = N,
                              n_chain = n_chain,
                              n_burn = n_burn,
                              n_samp = n_samp,
                              n_thin = n_thin,
                              prior = NULL,
                              DIC = TRUE,
                              parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}


#### Radiation model ####
fit.radiation.model <- function(M,s_ij, N, n_chain, n_burn, n_samp, n_thin, prior = NULL, DIC = TRUE, parallel = TRUE){
  
  fit.rad.model <-function (M, 
                            s_ij,
                            m_tot,
                            N = NULL, 
                            N_orig = NULL, 
                            N_dest = NULL, 
                            n_chain = 2, n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                            parallel = FALSE)
  {
    if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig <- N
    }
    else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
      N_dest <- N_orig
    }
    if (!(identical(dim(M)[1], dim(s_ij)[1], length(N_orig))))
      stop("Dimensions of input data must match")
    if (!(identical(dim(M)[2], dim(s_ij)[2], length(N_dest))))
      stop("Dimensions of input data must match")
    message(paste("::Fitting radiation model for", dim(M)[1],
                  "origins and", dim(M)[2], "destinations::",
                  sep = " "))
    if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
      M[, ] <- as.integer(M)
      N_orig[] <- as.integer(N_orig)
      N_dest[] <- as.integer(N_dest)
    }
    diag(M) <- 0
    if (is.null(prior)) {
      message("Using uniformative priors")
      null_prior <- c(1, 0.5)
      prior <- list(theta = null_prior)
    }
    else {
      message("Using supplied informative priors")
    }
    jags_data <- list(M = M, 
                      s_ij = s_ij,
                      m_tot = m_tot,
                      N_orig = N_orig, 
                      N_dest = N_dest,
                      prior_theta = prior$theta)
    
    jags_model <- "  
  model {
  
      # Poisson likelihood
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
  
          M[i,j] ~ dpois(c[i,j])
        }
  
      }
  
      # Radiation model
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
  
          c[i,j] <- ifelse(
            i == j,
            1e-06,
           (N_orig[i] * theta / (1 - N_orig[i] / m_tot)) * (N_orig[i] * N_dest[j]) / ((N_orig[i] + s_ij[i,j])*(N_orig[i] + N_dest[j] + s_ij[i,j]))
  
          )
        }
      }
  
      # Priors
      theta ~ dgamma(prior_theta[1], prior_theta[2]) # <- 1
  
    }"
    
    params <- c("theta")
    
    fit_jags(jags_data = jags_data, jags_model = jags_model,
             params = params, n_chain = n_chain, n_burn = n_burn,
             n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
  }
  output <- fit.rad.model(M = monthly.trips,
                          s_ij,
                          m_tot, 
                          N = N,
                          n_chain = n_chain,
                          n_burn = n_burn,
                          n_samp = n_samp,
                          n_thin = n_thin,
                          prior = NULL,
                          DIC = TRUE,
                          parallel = TRUE)
  
  parameters <- mobility::summarize_mobility(output)
  return(parameters)
}



#### Compiling Trip Estimates ####

trip.metrics <- function(x, M.predicted, M.observed, model, district.label, district.ID){
  
  row.names(x) <- colnames(x) <- as.integer(district.ID)
  M.predicted <- reshape2::melt(x, 
                                id.vars = c("origin", "destination"))  # prepare matrices for joining for chi square test
  M.predicted$model <- rep(model, nrows = length(M.predicted))
  colnames(M.predicted) <- c("origin", "destination","Trips.pred", "model")
  
  M.observed.long <- reshape2::melt(M.observed, 
                                    id.vars = c("origin", "destination"))  # prepare matrices for joining for chi square test
  M.observed.long$model <- rep("Observed", nrows = length(M.observed.long))
  colnames(M.observed.long) <- c("origin", "destination","Trips.obs", "model")
  
  M.compare <- left_join(M.predicted, M.observed.long, by = c("origin", "destination"))
  M.compare <- left_join(M.compare, district.label, by = c('origin'='adm2.code'))
  M.compare <- left_join(M.compare, district.label, by = c('destination'='adm2.code'))
  colnames(M.compare) <- c('origin', 'destination', 'Trips.pred', 'model.x', 'Trips.obs', 'model.y', 'origin.idx', 'destination.idx')
  
  M.compare$over.under.est.ratio <- (M.compare$Trips.pred + 0.1)/(M.compare$Trips.obs + 0.1) # added in 0.1 to avoid division by 0
  
  return(M.compare)
}


# predict trip numbers from models 1 and 2
sim.basic.gravity <- function (dist.fx, M.observed, N, D, theta = 1, omega.1 = 1, omega.2 = 1, gam = 1, model, district.ID, district.label = rank.Id, counts = FALSE) {
  if (!(identical(length(N), dim(D)[1], dim(D)[2])))
    stop("Check dimensions of input data N and D")
  if (!(length(c(theta, omega.1, omega.2)) == 3))
    stop("theta and omega parameters must be scalars")
  n.districts <- length(N)
  x <- f.d <- matrix(NA, n.districts, n.districts)
  error.fit.origin <- matrix(NA, n.districts, 4)
  
  for (i in 1:n.districts) {
    for (j in 1:n.districts) {
      if (i == j) {
        x[i, j] <- 0
      }
      else {
        if (length(gam) == 1 & dist.fx == "pwr") {
          f.d[i, j] <- D[i,j]^gam
        }
        else if (length(gam) == 1 & dist.fx == "exp") {
          f.d[i, j] <- exp(D[i,j]/gam)
        }
        else {
          stop("Incorrect length for parameters in Gamma dispersal kernel(s)")
        }
        x[i, j] <- exp(log(theta) + (omega.1 * log(N[i]) +
                                       omega.2 * log(N[j]) - log(f.d[i, j])))
      }
    }
  }
  
  M.compare <- trip.metrics(x, M.predicted, M.observed, model, district.label, district.ID)
  return(M.compare)
}


# predict trip numbers from models tripType
sim.trip.type.gravity <- function (dist.fx, M.observed, N, D, TT, theta = 1, omega.1 = 1, omega.2 = 1, gam = 1, model, district.ID, district.label = rank.Id, counts = FALSE){
  if (!(identical(length(N), dim(D)[1], dim(D)[2])))
    stop("Check dimensions of input data N and D")
  n.districts <- length(N)
  x <- f.d <- matrix(NA, n.districts, n.districts)
  error.fit.origin <- matrix(NA, n.districts, 4)
  
  for (i in 1:n.districts) {
    for (j in 1:n.districts) {
      n <- TT[i,j]
      if (i == j) {
        x[i, j] <- 0
      }
      else {
        if (dist.fx == "pwr") {
          f.d[i, j] <- D[i,j]^gam[n]
        }
        else if (dist.fx == "exp") {
          f.d[i, j] <- exp(D[i,j]/gam[n])
          
        }
        else {
          stop("Incorrect length for parameters in Gamma dispersal kernel(s)")
        }
        x[i, j] <- exp(log(theta) + (omega.1[n] * log(N[i]) +
                                       omega.2[n] * log(N[j]) - log(f.d[i, j])))
      }
    }
  }
  M.compare <- trip.metrics(x, M.predicted, M.observed, model, district.label, district.ID)
  return(M.compare)
}

# predict trip numbers from radiation model
sim.rad.mod <- function (M.observed, N, s_ij, theta, M, model, district.ID, district.label = rank.Id, counts = FALSE){
  if (!(identical(length(N), dim(s_ij)[1], dim(s_ij)[2])))
    stop("Check dimensions of input data N and s.ij")
  n.districts <- length(N)
  x <- matrix(NA, n.districts, n.districts)
  error.fit.origin <- matrix(NA, n.districts, 4)
  
  for (i in 1:n.districts) {
    
    for (j in 1:n.districts) {
      
      x[i,j] <- ifelse(
        i == j,
        0,
        (N[i] * theta / (1 - N[i] / M)) * (N[i] * N[j]) / ((N[i] + s_ij[i,j])*(N[i] + N[j] + s_ij[i,j]))
      )
    }
  }
  
  M.compare <- trip.metrics(x, M.predicted, M.observed, model, district.label, district.ID)
  return(M.compare)
}

# predict trip numbers from radiation model
sim.rad.mod <- function (M.observed, N, s_ij, theta, m_tot, model, district.ID, district.label = rank.Id, counts = FALSE){
  if (!(identical(length(N), dim(s_ij)[1], dim(s_ij)[2])))
    stop("Check dimensions of input data N and s.ij")
  n.districts <- length(N)
  x <- matrix(NA, n.districts, n.districts)
  error.fit.origin <- matrix(NA, n.districts, 4)
  
  for (i in 1:n.districts) {
    for (j in 1:n.districts) {
      
      x[i,j] <- ifelse(
        i == j,
        0,
        (N[i] * theta / (1 - N[i] / m_tot)) * (N[i] * N[j]) / ((N[i] + s_ij[i,j])*(N[i] + N[j] + s_ij[i,j]))
      )
    }
  }
  
  M.compare <- trip.metrics(x, M.predicted, M.observed, model, district.label, district.ID)
  return(M.compare)
}