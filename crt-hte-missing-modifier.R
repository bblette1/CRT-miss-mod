# Exploratory simulation of HTE in CRT w/ missing moderator
# Load libraries
library(dplyr)
library(geepack)
library(ggplot2)

# Simulator function
simulator <- function(trial, ICC_out, ICC_mod, ICC_miss, num_clusters) {
  
  # Data simulation
  size_clusters <- rpois(num_clusters, 50)
  
  # Create data frame
  df <- data.frame(cluster_ID = rep(1:num_clusters, size_clusters))
  df <- df %>% group_by(cluster_ID) %>%
    mutate(ind_ID = row_number(cluster_ID))
  
  # Treatment arm
  treated <- sample(1:num_clusters, num_clusters/2, replace = FALSE)
  df$A <- 1*(df$cluster_ID %in% treated)
  
  # Covariate
  n <- dim(df)[1]
  df$X <- rnorm(n, 0, 1)
  
  # Modifier
  alpha0_var <- pi^2 * ICC_mod / 3 / (1 - ICC_mod)
  alpha0 <- rnorm(num_clusters, 0, sqrt(alpha0_var))
  Mlogit <- 0.5 + alpha0
  Mprob <- exp(Mlogit) / (1 + exp(Mlogit))
  Mfull <- rbinom(n, 1, Mprob)
  
  # Outcome
  Y_resvar <- 3
  alpha1_var <- Y_resvar * ICC_out / (1 - ICC_out)
  alpha1 <- rnorm(num_clusters, 0, sqrt(alpha1_var))
  df$Y <- 1 + 1.5*df$A + df$Mfull - 0.75*df$A*df$Mfull + 0.5*df$X +
    rep(alpha1, size_clusters) + rnorm(n, 0, Y_resvar)
  
  # Missingness
  Rlogit <- 0.5 + 0.5*df$X
  Rprob <- exp(Rlogit) / (1 + exp(Rlogit))
  df$R <- rbinom(n, 1, Rprob)
  
  df$M <- NA
  df$M[df$R == 1] <- df$Mfull[df$R == 1]
  
}