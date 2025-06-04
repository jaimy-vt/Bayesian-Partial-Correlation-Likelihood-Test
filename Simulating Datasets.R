library(Matrix)
library(igraph)
library(corpcor)
library(bootnet)
library(qgraph)
library(mvtnorm)
library(MASS)
library(lavaan)

#### Configurations
nobs <- seq(100, 2000, 100) # Different sample sizes per dataset
ndif <- 100
nds <- length(nobs)         # Number of simulated network datasets
nvs <- c(5, 10, 15)         # Number of variables
edge_dens <- c(0.8)         # Edge density for the networks

save_location_factor = '/factor_datasets'

#### Simulating data from a factor model ####

f_datasets <- c()

print("Starting factor dataset generating.")


for (v in 1:length(nvs)){
  nv <- nvs[v]
  f_datasets[[paste0("nv_", nv)]] <- list()
  
  for (i in 1:nds){
    n <- nobs[i]     # Number of observations
    f_datasets[[paste0("nv_", nv)]][[paste0("nobs_", n)]] <- list()
    
    for (di in 1:ndif){
      print(paste('Variables:', nvs[v], 'Samples size:', n, 'dataset', di, 'out of', ndif))
      
      # Step 1: Define factor loadings (these should typically be between 0 and 1)
      loadings <- runif(nv, 0.1, 0.9)
      signs <- sample(c(-1, 1), nv, replace = TRUE)  # Randomly assign -1 or 1
      loadings <- loadings * signs  # Apply the sign
      
      # Step 2: Simulate the latent factor (F)
      eta <- rnorm(n, mean = 0, sd = 1)  # Standard normal latent factor
      
      # Step 3: Simulate unique errors (assumed to be normally distributed)
      errorcov<-diag(0.5,nv)
      errors<-mvrnorm(n,mu=rep(0,nv),Sigma=errorcov)
      #errors <- matrix(rnorm(n * nv, mean = 0, sd = 0.5), nrow = n, ncol = nv)
      
      # Step 4: Generate observed variables (X1 to X5)
      X <- eta %*% t(loadings) + errors  # Matrix multiplication for factor structure
      
      # Convert to a data frame for easier handling
      factor_data <- as.data.frame(X)
      colnames(factor_data) <- paste0("V", 1:nv)  # Name the variables
      
      f_datasets[[paste0("nv_", nv)]][[paste0("nobs_", n)]][[di]] <- factor_data
    }
  }
}

##### Saving the simulated datasets for later use ####
save(f_datasets, file = save_location_factor)

#### Simulating Data from a Network ####
save_location_network = '/network_datasets'

n_datasets <- c() # List where the network simulated datasets will be stored

print("Starting network dataset generating.")
for (d in 1:length(edge_dens)){
  dens <- edge_dens[d]
  n_datasets[[paste0("dens_", dens)]] <- list()
  
  for (v in 1:length(nvs)){
    nv <- nvs[v]
    n_datasets[[paste0("dens_", dens)]][[paste0("nv_", nv)]] <- list()
    
    for (i in 1:nds){
      n_size <- nobs[i]
      n_datasets[[paste0("dens_", dens)]][[paste0("nv_", nv)]][[paste0("nobs_", n_size)]] <- list()
      
      for (di in 1:ndif){
        print(paste('Density:', edge_dens[d], 'Variables:', nvs[v], 'Samples size:', n_size, 'dataset', di, 'out of', ndif))
        
        repeat { # Making sure I only use the positive definite matrices
          
          # Keeping the bounds negative, so the connections are positive
          lower_bound <- -4
          upper_bound <- 4
          
          # Determining which off diagonal cells are not zero
          edge_density_matrix <- erdos.renyi.game(nv, dens) # This keeps the diagonal zero
          edm <- as.matrix(get.adjacency(edge_density_matrix))
          edm
          
          # Assigning a random number between -1 and 1 for the non zero cells
          precision <- matrix(runif(nv^2, lower_bound, upper_bound), nv, nv) * edm
          precision
          
          diag(precision) <- 0 # Making sure the diagonal stays 0
          precision
          
          # Making all the off diagonals a number between 0 and 1 (absolute)
          precision <- precision * 0.1
          precision
          
          
          # Changing the diagonal to the smallest eigenvalue + some arbitrary number
          diag(precision) <- abs(min(Re(eigen(precision)$values))) + 0.2
          precision
          
          # Force the matrix to be symmetric
          precision <- as.matrix(forceSymmetric(precision))
          precision
          is.positive.definite(precision)
          
          # Checking whether the matarix is positive definite
          if (is.positive.definite(precision)) {
            break
          }
        }
        
        # First, calculating the diagonal matrix with the formula delta_jj = p_jj ^-(1/2)
        diagonal <- diag(1/sqrt(diag(precision)))
        diagonal
        
        # Pre- and post- multiplying the precision matrix by the diagonal matrix
        standardised_precision <- diagonal %*% precision %*% diagonal
        standardised_precision
        
        # Subtracting this standardized precision matrix from an identity matrix
        identity <- diag(nv)
        weights <- identity - standardised_precision
        weights
        
        weights <- as.matrix(weights)
        diag(weights) <- 1 
        weights
        is.positive.definite(weights)
        
        # Calculating the correlation matric from the weights (partial correlations) matrix
        sigma <- pcor2cor(weights)
        is.positive.definite(sigma)
        
        #Simulating Data from the Weight Matrix
        sim_data <- mvrnorm(n_size, mu = rep(0, nv), Sigma = sigma)
        
        network_data <- as.data.frame(sim_data)
        colnames(network_data) <- paste0("V", 1:nv)  # Name the variables
        
        n_datasets[[paste0("dens_", dens)]][[paste0("nv_", nv)]][[paste0("nobs_", n_size)]][[di]] <- network_data
      }
    }
  }
}


##### Saving the simulated datasets for later use ####
title <- paste("dens", edge_dens, "network_datasets", sep = "_")
file_loc <- paste(save_location_network, title, sep = '/')
save(n_datasets, file = file_loc)

#### Simulating Data from Network and Factor Combined ####
"Creating a combined dataset by simulating data from a factor model that has
correlated errors. The correlated errors represents the network"

save_location_combined = '/combined_datasets'

nds <- 100
nv <- 10
Nobs <- seq(100, 2000, 100)
edge_dens <- 0.5
F_infl <- c(.2, .5, .8) # Whether factor of network has more influence on the data

for (finf in 1:length(F_infl)){
  f_infl <- F_infl[[finf]]
  combined_datasets <- c()
  
  for (obs in 1:length(Nobs)){
    nobs <- Nobs[obs]
    combined_datasets[[paste0("nobs_", nobs)]] <- list()
    for (i in 1:nds){
      n <- nobs
      loadings <- runif(nv, 0.1, 0.9)
      signs <- sample(c(-1, 1), nv, replace = TRUE)  # Randomly assign -1 or 1
      loadings <- loadings * signs  # Apply the sign
      
      errorcov <- diag((1 - (loadings^2)), nv)
      
      cor_matrix_factor <- loadings %*% t(loadings) + errorcov
      is.positive.definite(cor_matrix_factor)
      
      repeat{
        graph <- erdos.renyi.game(nv, edge_dens, directed = FALSE)
        adj_matrix <- as.matrix(get.adjacency(graph))
        
        lower_bound <- -4
        upper_bound <- 4
        
        net_pcor_matrix <- matrix(runif(nv^2, lower_bound, upper_bound), nv, nv) * adj_matrix
        diag(net_pcor_matrix) <- 0
        net_pcor_matrix <- net_pcor_matrix * 0.1
        diag(net_pcor_matrix) <- abs(min(Re(eigen(net_pcor_matrix)$values))) + 0.2
        net_pcor_matrix <- as.matrix(forceSymmetric(net_pcor_matrix))
        
        is.positive.definite(net_pcor_matrix)
        
        net_cor_matrix <- corpcor::pcor2cor(net_pcor_matrix)
        net_cor_matrix
        
        if (any(is.nan(net_cor_matrix))){
          warning("NaNs found in network correlation matrix. Restarting.")
          next
        } else if (!any(is.nan(net_cor_matrix))){
          cov_matrix <- f_infl * cor_matrix_factor + (1-f_infl) * net_cor_matrix
          cov_matrix
          #cor_matrix <- cov2cor(cov_matrix)
          #cor_matrix
          
          if (!is.positive.definite(cov_matrix)){
            warning("Resulting covariance matrix not positive definite. Restarting.")
            next
          } else if(is.positive.definite(cov_matrix)){
            break
          }
        }
      }
      
      sim_data <- mvrnorm(n = n, mu = rep(0, ncol(cov_matrix)), Sigma = cov_matrix)
      sim_data <- as.data.frame(sim_data)
      head(sim_data)
      
      combined_datasets[[paste0("nobs_", nobs)]][[i]] <- sim_data
    }
  }
  
  title <- paste('f_infl', f_infl, 'combined_datasets', sep = "_")
  file_loc <- paste(save_location_combined, title, sep = '/')
  save(combined_datasets, file = file_loc)
}



