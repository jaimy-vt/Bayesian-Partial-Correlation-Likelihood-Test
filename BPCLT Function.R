#### Loading in the packages ####
library(MASS)
library(blavaan)
library(coda)
library(easybgm)
library(corpcor)
library(lavaan)
library(parallel)
library(future.apply)

n_iter <- 1000

#### BPCLT Function ####
BPCLT <- function(dataset) {
  
  netlikely <- c()
  falikely <- c()
  
  n_ident <- 0
  
  MyData <- dataset
  nv <- ncol(MyData)
  
  "Factor fit"
  # For every dataset, calculate the posterior density of the model parameters here
  factor_model <- paste("procor =~", paste0("V", 1:nv, collapse = " + "))
  
  factor_fit <- bsem(factor_model, data = MyData, std.lv = TRUE)
  summary(factor_fit)
  
  est_factor_loadings <- blavInspect(factor_fit, what = 'mcmc')[,1:nv]
  combined_samples <- as.mcmc.list(est_factor_loadings)
  posterior_loadings <- do.call(rbind, combined_samples)
  
  est_variances <- blavInspect(factor_fit, what = 'mcmc')[,(nv+1):(nv+nv)]
  combined_variances <- as.mcmc.list(est_variances)
  posterior_variances <- do.call(rbind, combined_variances)
  
  "Network Fit"
  network_fit <- easybgm::easybgm(data = MyData, type = 'continuous', 
                                  save = TRUE, package = 'BDgraph')
  
  net_sum <- summary(network_fit)
  
  right_order <- order(net_sum$parameters$Relation)
  category <- net_sum$parameters$Category[right_order]
  
  posterior_weights <- network_fit$samples_posterior
  
  # If easybgm says excluded, set all values of posterior density to 0
  for (c in 1:length(category)){
    if(category[c] == 'excluded'){
      posterior_weights[,c] <- 0
    }
  }
  
  
  # Calculating the anomalous correlation proportion for this dataset
  "Anomalous Correlation in the dataset"
  cor_data <- cor(MyData)
  
  pcor_data <- cor2pcor(cor_data)
  U_data_cor <- cor_data[lower.tri(cor_data, diag = F)]
  
  #unique elements in correlation matrix
  U_data_pcor <- pcor_data[lower.tri(pcor_data, diag = F)] #unique elements in partial correlation matrix
  
  U_FA_pcor <- pcor_data[lower.tri(pcor_data, diag = F)]
  U_FA_cor <- cor_data[lower.tri(pcor_data, diag = F)]
  
  signswitch_data <- which(U_data_pcor * U_data_cor < 0) #number of correlations that have different sign than pcor
  stronger_data <- which((U_FA_pcor^2) >(U_FA_cor^2))
  together_data <- union(signswitch_data, stronger_data)
  total_data <- length(together_data)
  
  propData <- total_data/length(U_data_cor) #proportion of correlations that have different sign than pcor
  
  propFA <- rep(NA, n_iter)
  propNW <- rep(NA, n_iter)
  
  for (j in 1:n_iter){
    print(paste("Iteration", j, '/', n_iter))
    
    "Factor Dataset"
    # For every iteration, randomly sample from the model parameters density
    # distribution and use that model to generate a new dataset
    sampled_row_f <- sample(1:nrow(posterior_loadings), 1)
    loadings <- posterior_loadings[sampled_row_f,]
    variances <- posterior_variances[sampled_row_f,]
    
    
    # Simulating data from these factor loadings
    lambda <- matrix(loadings, length(loadings))
    theta <- diag(variances)
    
    sigma_factor <- lambda %*% t(lambda) + theta
    
    DataFA <- mvrnorm(nrow(MyData), mu = rep(0, nrow(sigma_factor)), Sigma = sigma_factor)
    
    "Network Dataset"
    attempt <- 0
    PDattempt <- 0
    nan_in_nearPD_cor <- 0
    repeat{
      attempt <- attempt + 1
      sampled_row_n <- sample(1:nrow(posterior_weights), 1)
      weights <- posterior_weights[sampled_row_n, ]
      
      # Simulating data from these network weights
      pcor_network <- matrix(0, nrow = nv, ncol = nv)
      index_upper_triangle <- which(lower.tri(pcor_network), arr.ind = TRUE)
      
      # Assign the weights to the corresponding upper triangle entries
      for (ix in 1:length(index_upper_triangle[,1])) {
        row <- index_upper_triangle[ix, 1]
        col <- index_upper_triangle[ix, 2]
        pcor_network[row, col] <- weights[ix]
      }
      
      pcor_network <- pcor_network + t(pcor_network)
      diag(pcor_network) <- 1
      
      if (is.positive.definite(pcor_network) & attempt < 10) {
        cor_network <- pcor2cor(pcor_network)
        if (any(is.nan(cor_network))) {
          next  # If NaNs are present, skip this iteration and retry
        } else {
          sigma_network <- cor_network
          if (is.positive.definite(sigma_network)){
            print("Positive Definite sigma. No nearPD attempt")
            break
          } else {
            print("Sigma not PD. Next attempt")
            next
          }
        }
      } else if(!is.positive.definite(pcor_network) & attempt < 10) {
        warning("Pcor is not PD. Skipping to next iteration")
        next
      }
      
      if (attempt >= 10) {
        warning(paste("Max attempts reached. Using nearPD to approximate covariance matrix. Attempt:", PDattempt, "/10"))
        PDattempt <- PDattempt + 1
        pcor_network <- as.matrix(Matrix::nearPD(pcor_network, corr = TRUE)$mat)
        cor_network <- pcor2cor(pcor_network)
        
        if (any(is.nan(cor_network)) & PDattempt < 10){
          warning("NaNs in nearPD cor Matrix. Skipping to next iteration")
          nan_in_nearPD_cor <- nan_in_nearPD_cor + 1
          next
        } else if (!any(is.nan(cor_network)) & PDattempt < 10){
          if(is.positive.definite(cor_network)){
            sigma_network <- cor_network
            break
          } else {
            next
          }
        }
        
        if (PDattempt >= 10){
          warning("nearPD still resulted in NaNs. Using identity matrix as fallback.")
          sigma_network <- diag(3)  # Fallback to identity matrix
          n_ident <- n_ident + 1
          break
        } else {
          next
        }
      }
    }
    
    DataNW <- mvrnorm(nrow(MyData), mu = rep(0, nrow(sigma_network)), Sigma = sigma_network)
    
    # Calculate the anomalous correlation proportion from this dataset and save
    corFA <- cor(DataFA)
    pcorFA <- cor2pcor(corFA)
    
    corNW <- cor(DataNW)
    pcorNW <- cor2pcor(corNW)
    
    U_FA_cor <- corFA[lower.tri(corFA, diag = F)]
    U_FA_pcor <- pcorFA[lower.tri(pcorFA, diag = F)]
    
    U_NW_cor <- corNW[lower.tri(corNW, diag = F)]
    U_NW_pcor <- pcorNW[lower.tri(pcorNW, diag = F)]
    
    signswitch_FA <- which(U_FA_pcor * U_FA_cor < 0)
    signswitch_NW <- which(U_NW_pcor * U_NW_cor < 0)
    
    stronger_FA <- which((U_FA_pcor ^ 2) >(U_FA_cor^2))
    stronger_NW <- which((U_NW_pcor ^ 2) > (U_NW_cor^2))
    
    together_FA <- union(signswitch_FA, stronger_FA)
    together_NW <- union(signswitch_NW, stronger_NW)
    
    total_FA <- length(together_FA)
    total_NW <- length(together_NW)
    
    propFA[j] <- total_FA/length(U_FA_cor)
    propNW[j] <- total_NW/length(U_NW_cor)
  }
  
  # Compute density objects
  density_NW <- density(propNW, bw = 0.1)
  density_FA <- density(propFA, bw = 0.1)
  
  density_net <- approxfun(density_NW$x, density_NW$y, rule = 2)
  density_fac <- approxfun(density_FA$x, density_FA$y, rule = 2)
  
  point <- propData
  
  netdens <- density_net(point)
  netlikely <- netdens
  fadens <- density_fac(point)
  falikely <- fadens
  
  results <- list(network_dens = density_NW,
                  factor_dens = density_FA,
                  identity = n_ident)
  
  return(results)
}

#### future.apply parralel ####

# Check the number of available cores
available_cores <- parallel::detectCores()
print(paste("Available cores:", available_cores))

# Use more workers (but leave 1 core free for system processes)
plan(multisession, workers = available_cores - 1)
results <- future_lapply(datasets, BPCLT, future.seed = TRUE)

