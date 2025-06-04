
library(MASS)
library(blavaan)
library(coda)
library(easybgm)
library(corpcor)
library(lavaan)

#### First, loading in the pre-generated datasets ####
load(file = "dens_0.2_network_datasets")
network_dens_02 <- n_datasets
load(file = "dens_0.5_network_datasets")
network_dens_05 <- n_datasets
load(file = "dens_0.8_network_datasets")
network_dens_08 <-n_datasets

#### Organizing datasets into a list for easier selection ####
n_datasets <- list(
  "0.2" = network_dens_02$dens_0.2,
  "0.5" = network_dens_05$dens_0.5,
  "0.8" = network_dens_08$dens_0.8
)

# Loading in the factor datasets
load(file = "/factor_datasets")

#### Specifying which dataset characteristics to analyze #### 
true_model <- "factor" # "factor" or "network"
nobs <- "nobs_100"
nvs <- "nv_5"
edge_dens <- "0.5"
n_iter <- 1000

#### PCLT on the selected dataset ####
start_time <- Sys.time()

netlikely <- c()
falikely <- c()

n_ident <- 0

##### Network Datasets Analysis #####

if (true_model == "network"){
  network_datasets <- n_datasets[[edge_dens]][[nvs]][[nobs]][1:5]
  length(network_datasets)
  
  n_correct <- 0
  n_incorrect <- 0
  
  for (i in 1:length(network_datasets)){
    MyData <- network_datasets[[i]]
    nv <- ncol(MyData)
    
    print(paste("Dataset", i, "out of", length(network_datasets)))
    
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
      print(paste("Dataset", i, '/', length(network_datasets), "Iteration", j, '/', n_iter))
      
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
    netlikely[i] <- netdens
    fadens <- density_fac(point)
    falikely[i] <- fadens
    
    if (netdens > fadens){
      n_correct <- n_correct + 1
    } else if (fadens > netdens){
      n_incorrect <- n_incorrect + 1
    }
  }
}

##### Factor Datasets Analysis #####
if (true_model == "factor"){
  factor_datasets <- f_datasets[[nvs]][[nobs]][1:5]
  
  f_correct <- 0
  f_incorrect <- 0
  
  for (i in 1:length(factor_datasets)){
    MyData <- factor_datasets[[i]]
    nv <- ncol(MyData)
    
    print(paste("Dataset", i, "out of", length(factor_datasets)))
    
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
      print(paste("Dataset", i, '/', length(factor_datasets), "Iteration", j, '/', n_iter))
      
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
    netlikely[i] <- netdens
    fadens <- density_fac(point)
    falikely[i] <- fadens
    
    if (fadens > netdens){
      f_correct <- f_correct + 1
    } else if (netdens > fadens){
      f_incorrect <- f_incorrect + 1
    }
  }
}

end_time <- Sys.time()

#### Proportion correct ####
time_taken <- end_time - start_time
time_taken

if (true_model == "factor"){
  prop_correct <- f_correct/length(factor_datasets)
  print(paste("Proportion correct:", prop_correct))
  prop_incorrect <- f_incorrect/length(factor_datasets)
  print(paste("Proportion incorrect:", prop_incorrect))
}

if (true_model == "network"){
  prop_correct <- n_correct/length(network_datasets)
  print(paste("Proportion correct:", prop_correct))
  prop_incorrect <- n_incorrect/length(network_datasets)
  print(paste("Proportion incorrect:", prop_incorrect))
}

print(paste("Calculation time:", time_taken))

# Create a named vector with the proportions
proportions <- c(Correct = prop_correct, Incorrect = prop_incorrect)

if (true_model == "network"){
  xlb <- paste("Observed model:", true_model, edge_dens)
} else if (true_model == "factor"){
  xlb <- paste("Observed model:", true_model)
}

# Create the bar plot and store the bar midpoints
bar_midpoints <- barplot(proportions, 
                         main = "Proportion of Correct vs Incorrect",
                         col = c("forestgreen", "darkred"), 
                         ylim = c(0, 1),  # Adjust as needed
                         ylab = "Proportion",
                         xlab = xlb)

# Add numbers on top of the bars
text(x = bar_midpoints, y = proportions + 0.1, 
     labels = round(proportions, 2), cex = 1.2, col = "black", font = 2)
text(x = 2.25, y = .9, 
     labels = paste(nobs, nvs), 
     adj = 1, cex = 0.7, col = "black", font = 2)


#### Plotting the last dataset ####
xlim_range <- range(c(density_NW$x, density_FA$x))
ylim_range <- range(c(density_NW$y, density_FA$y))

if (true_model == "factor"){
  title <- 'Last dataset example plot - Observed = factor'
} else if (true_model == "network"){
  title <- 'Last dataset example plot - Observed = network'
}

plot(density_NW, col = 'blue', main = title,
     xlim = xlim_range, ylim = ylim_range)
lines(density_FA, col = 'red')
abline(v = point, lty = "dashed")
legend("topright", legend = c('Network', 'Factor', 'Observed'), 
       col = c('blue', 'red', 'black'), lty = c('solid', 'solid', 'dashed'))

if (true_model == 'network'){
  print(table(netlikely > falikely))
  print(mean(netlikely - falikely))
} else if (true_model == 'factor'){
  print(table(falikely > netlikely))
  print(mean(falikely - netlikely))
}

#### Warnings ####
summary(warnings())

print(paste("Number of identity matrices used for network data simulation:", n_ident))


#### Saving the results into a file ####
if (true_model == 'network'){
  filename <- paste('network', 'dens', edge_dens, nvs, nobs, sep = '_')
} else if (true_model == 'factor'){
  filename <- paste('factor', nvs, nobs, sep = '_')
}

save_file <- paste("~/UvA/Master/Scriptie/Code/Bayesian Test/BPCLT/Seed Results/", filename, sep = '')
#save.image(file = save_file)
#load(save_file)

#### BPCLT Functions ####

"The only difference between these two functions is how correct and incorrect is specified "
BPCLT_network <- function(dataset) {
  n_correct <- 0
  n_incorrect <- 0
  
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
  
  if (netdens > fadens){
    n_correct <- n_correct + 1
  } else if (fadens > netdens){
    n_incorrect <- n_incorrect + 1
  }
  
  results <- list(correct = n_correct, 
                  incorrect = n_incorrect,
                  network_dens = density_NW,
                  factor_dens = density_FA,
                  identity = n_ident)
  
  return(results)
}

BPCLT_factor <- function(dataset) {
  f_correct <- 0
  f_incorrect <- 0
  
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
  
  if (fadens > netdens){
    f_correct <- f_correct + 1
  } else if (netdens > fadens){
    f_incorrect <- f_incorrect + 1
  }
  
  results <- list(correct = f_correct, 
                  incorrect = f_incorrect,
                  network_dens = density_NW,
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

if (true_model == 'network'){
  datasets <- n_datasets[[edge_dens]][[nvs]][[nobs]]
  start_time_par <- Sys.time()
  results <- future_lapply(datasets, BPCLT_network, future.seed = TRUE)
  end_time_par <- Sys.time()
} else if (true_model == 'factor'){
  datasets <- f_datasets[[nvs]][[nobs]]
  start_time_par <- Sys.time()
  results <- future_lapply(datasets, BPCLT_factor, future.seed = TRUE)
  end_time_par <- Sys.time()
}

# Aggregate results
sum(sapply(results, function(res) res$identity))

if (true_model == 'network'){
  n_correct_total <- sum(sapply(results, function(res) res$correct))
  n_incorrect_total <- sum(sapply(results, function(res) res$incorrect))
  
  print(paste("Total network correct:", n_correct_total))
  print(paste("Total network incorrect:", n_incorrect_total))
} else if (true_model == 'factor'){
  f_correct_total <- sum(sapply(results, function(res) res$correct))
  f_incorrect_total <- sum(sapply(results, function(res) res$incorrect))
  
  print(paste("Total factor correct:", f_correct_total))
  print(paste("Total factor incorrect:", f_incorrect_total))
}

#### Plotting the proportions correct ####

if (true_model == "network"){
  proportions <- c(Correct = n_correct_total/100, Incorrect = n_incorrect_total/100)
} else if (true_model == 'factor'){
  proportions <- c(Correct = f_correct_total/100, Incorrect = f_incorrect_total/100)
}

if (true_model == "network"){
  xlb <- paste("Observed model:", true_model, edge_dens)
} else if (true_model == "factor"){
  xlb <- paste("Observed model:", true_model)
}

# Create the bar plot and store the bar midpoints
bar_midpoints <- barplot(proportions, 
                         main = "Proportion of Correct vs Incorrect",
                         col = c("forestgreen", "darkred"), 
                         ylim = c(0, 1.1),  # Adjust as needed
                         ylab = "Proportion",
                         xlab = xlb)

# Add numbers on top of the bars
text(x = bar_midpoints, y = proportions + 0.05, 
     labels = proportions, cex = 1, col = "black", font = 2)
text(x = 2.25, y = .9, 
     labels = paste(nobs, nvs), 
     adj = 1, cex = 0.7, col = "black", font = 2)


#### Saving the results into a file ####
if (true_model == 'network'){
  filename <- paste('network', 'dens', edge_dens, nvs, nobs, sep = '_')
} else if (true_model == 'factor'){
  filename <- paste('factor', nvs, nobs, sep = '_')
}

save_file <- paste(save_location, filename, sep = '')
#save.image(file = save_file)

