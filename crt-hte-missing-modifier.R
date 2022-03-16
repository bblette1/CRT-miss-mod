# Exploratory simulation of HTE in CRT w/ missing moderator
rm(list = ls())

# Load libraries
library(BayesLogit)
library(dbarts)
library(dplyr)
library(geepack)
library(ggplot2)
library(jomo)
library(lme4)
library(matrixStats)
library(MCMCpack)
library(mitml)
library(mvtnorm)

# Simulator function
simulator <- function(trial, ICC_out, ICC_mod, num_clusters) {
  
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
  Mlogit <- 0.5 + rep(alpha0, size_clusters)
  Mprob <- exp(Mlogit) / (1 + exp(Mlogit))
  df$Mfull <- rbinom(n, 1, Mprob)
  
  # Outcome
  Y_resvar <- 3
  alpha1_var <- Y_resvar * ICC_out / (1 - ICC_out)
  alpha1 <- rnorm(num_clusters, 0, sqrt(alpha1_var))
  df$Y <- 1 + 1.5*df$A + df$Mfull - 0.75*df$A*df$Mfull +
    #0.7*df$X +
    #0.7*df$X*df$A +
    #0.7*df$X*df$Mfull +
    #0.7*df$X*df$Mfull*df$A +
    0.8*df$X*df$A - 0.4*df$X*df$Mfull + 0.7*df$X*df$Mfull*df$A +
    rep(alpha1, size_clusters) + rnorm(n, 0, sqrt(Y_resvar))
  
  # Missingness
  Rlogit <- 1.2 + 0.5*df$X - 0.2*df$Y
  Rprob <- exp(Rlogit) / (1 + exp(Rlogit))
  df$R <- rbinom(n, 1, Rprob)
  
  df$M <- NA
  df$M[df$R == 1] <- df$Mfull[df$R == 1]
  
  # Add interactions to data for JAV analyses
  df$AM <- df$A*df$M
  df$XAM <- df$X*df$A*df$M
  
  
  # Data analysis
  # Truth
  truemod <- geeglm(Y ~ A*Mfull + X:A:Mfull, family = "gaussian", data = df,
                    id = cluster_ID, corstr = "exchangeable")
  
  # Method 1: CRA-GEE
  mod1 <- geeglm(Y ~ A*M, family = "gaussian",
                 data = df[!is.na(df$M), ],
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 2: Single Imputation - Passive
  impmod2_M <- glm(M ~ X*A*Y, data = df[!is.na(df$M), ],
                   family = "binomial")
  dfimp2 <- df
  dfimp2$M[is.na(dfimp2$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod2_M, df[is.na(df$M), ],
                                              type = "response"))
  mod2 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp2,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 3: Single Imputation - JAV
  impmod3_AM <- glm(AM ~ X*A*Y, data = df[!is.na(df$M), ],
                    family = "binomial")
  impmod3_XAM <- glm(XAM ~ X*A*Y, data = df[!is.na(df$M), ],
                     family = "gaussian")
  
  dfimp3 <- df
  dfimp3$M[is.na(dfimp3$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod2_M, df[is.na(df$M), ],
                                              type = "response"))
  dfimp3$AM[is.na(dfimp3$AM)] <- rbinom(sum(is.na(df$AM)), 1,
                                        predict(impmod3_AM, df[is.na(df$AM), ],
                                                type = "response"))
  #dfimp3$XAM[is.na(dfimp3$XAM)] <- predict(impmod3_XAM, df[is.na(df$XAM), ])
  mod3 <- geeglm(Y ~ A + M + AM, family = "gaussian", data = dfimp3,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 4: Single Imputation - JAV enhanced
  impmod4_AM <- glm(AM ~ X*Y, data = df[!is.na(df$M) & df$A != 0, ],
                    family = "binomial")
  
  dfimp4 <- df
  dfimp4$M[is.na(dfimp4$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod2_M, df[is.na(df$M), ],
                                              type = "response"))
  dfimp4$AM[is.na(dfimp4$AM) & dfimp4$A == 0] <- 0
  dfimp4$AM[is.na(dfimp4$AM) & dfimp4$A != 0] <-
    rbinom(sum(is.na(df$AM) & df$A != 0), 1,
           predict(impmod4_AM, df[is.na(df$AM) & df$A != 0, ],
                   type = "response"))
  mod4 <- geeglm(Y ~ A + M + AM, family = "gaussian", data = dfimp4,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Methods 3-8: Multiple Imputation methods
  numimp <- 15
  ests3 <- ests4 <- ests5 <- ests6 <- ests7 <- ests8 <-
    ests9 <- ests10 <- array(NA, dim = c(numimp, 4))
  varests3 <- varests4 <- varests5 <- varests6 <- varests7 <-
    varests8 <- varests9 <- varests10 <- array(NA, dim = c(numimp, 4))
  
  impmod4 <- glmer(M ~ X + A + Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  #impmod3 <- geeglm(M ~ X + A*Y, family = binomial(),
  #data = df[!is.na(df$M), ],
  #id = cluster_ID, corstr = "exchangeable")
  impmod5 <- glmer(M ~ X + A*Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  impmod6 <- glmer(M ~ X*A + Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  impmod7 <- glmer(M ~ X*A + Y*A + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  impmod8 <- glmer(M ~ X*A*Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  
  #impmodbart <- gbart(x.train = as.matrix(df[!is.na(df$M), c(3, 4, 6)]),
                      #y.train = as.vector(df$M[!is.na(df$M)]),
                      #type = "pbart")
  #predtest2 <- predict(impmodbart, as.matrix(df[is.na(df$M),c(3, 4, 6)]))
  
  #impmodbart <- rbart_vi(M ~ . - cluster_ID, n.trees = 75,
                         #df[!is.na(df$M), c(1, 3, 4, 6, 8)],
                         #group.by = cluster_ID, keepTrees = TRUE,
                         #n.chains = 2, power = 2)
  #predtest <- predict(impmodbart,
                      #df[is.na(df$M), c(1, 3, 4, 6, 8)],
                      #group.by = df[is.na(df$M),
                                    #c(1, 3, 4, 6, 8)]$cluster_ID,
                      #type = "ppd")
  
  for (m in 1:numimp) {
    
    # Method 3: Multiple Imputation
    dfimp3 <- df
    dfimp3$M[is.na(dfimp3$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod2_M, df[is.na(df$M), ],
             type = "response"))
    mod3 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp3,
                   id = cluster_ID, corstr = "exchangeable")
    ests3[m, ] <- coef(mod3)
    varests3[m, ] <- (summary(mod3)$coefficients[, 2])^2
    
    # Method 4: Multilevel MI cond. on X
    dfimp4 <- df
    dfimp4$M[is.na(dfimp4$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod4, df[is.na(df$M), ],
                                          type = "response"))
    mod4 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp4,
                   id = cluster_ID, corstr = "exchangeable")
    ests4[m, ] <- coef(mod4)
    varests4[m, ] <- (summary(mod4)$coefficients[, 2])^2
    
    # Method 5: Multilevel MI cond. on X and Y
    dfimp5 <- df
    dfimp5$M[is.na(dfimp5$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod5, df[is.na(df$M), ],
                                          type = "response"))
    mod5 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp5,
                   id = cluster_ID, corstr = "exchangeable")
    ests5[m, ] <- coef(mod5)
    varests5[m, ] <- (summary(mod5)$coefficients[, 2])^2
    
    # Method 6: Multilevel MI cond. on X, Y, and A
    dfimp6 <- df
    dfimp6$M[is.na(dfimp6$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod6, df[is.na(df$M), ],
                                          type = "response"))
    mod6 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp6,
                   id = cluster_ID, corstr = "exchangeable")
    ests6[m, ] <- coef(mod6)
    varests6[m, ] <- (summary(mod6)$coefficients[, 2])^2
    
    # Method 7: Multilevel MI cond. on X and A*Y
    dfimp7 <- df
    dfimp7$M[is.na(dfimp7$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod7, df[is.na(df$M), ],
                                          type = "response"))
    mod7 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp7,
                   id = cluster_ID, corstr = "exchangeable")
    ests7[m, ] <- coef(mod7)
    varests7[m, ] <- (summary(mod7)$coefficients[, 2])^2
    
    # Method 8: Multilevel MI cond. on X*Y*A
    dfimp8 <- df
    dfimp8$M[is.na(dfimp8$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod8, df[is.na(df$M), ],
                                          type = "response"))
    mod8 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp8,
                   id = cluster_ID, corstr = "exchangeable")
    ests8[m, ] <- coef(mod8)
    varests8[m, ] <- (summary(mod8)$coefficients[, 2])^2
    
    # Method 9: BART imputation
    #dfimp9 <- df
    #dfimp9$M[is.na(dfimp9$M)] <-
      #rbinom(sum(is.na(df$M)), 1, colMeans(predtest))
    #mod9 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp9,
                   #id = cluster_ID, corstr = "exchangeable")
    #ests9[m] <- coef(mod9)[4]
    #varests9[m] <- (summary(mod9)$coefficients[4, 2])^2
    
  }
  
  # Data augmentation
  #dfimp10 <- df
  #dfimp10$M[is.na(dfimp10$M)] <-
    #rbinom(sum(is.na(dfimp10$M)), 1, mean(dfimp10$M, na.rm = T))
  #ests10 <- array(NA, dim = c(100, 5))
  #varests10 <- array(NA, dim = c(100, 5))
  #for (g in 1:100) {
    
    #impmod10 <- glmer(M ~ X*A*Y + (1 | cluster_ID), family = "binomial",
                     #data = dfimp10)
    #mod10 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp10,
                    #id = cluster_ID, corstr = "exchangeable")
    #ests10[g, ] <- coef(mod10)
    #varests10[g, ] <- (summary(mod10)$coefficients[, 2])^2
    
    #dfimp10 <- df
    #dfimp10$M[is.na(dfimp10$M)] <-
      #rbinom(sum(is.na(df$M)), 1, predict(impmod10, df[is.na(df$M), ],
                                          #type = "response"))
    
  #}
  
  # Use jomo package
  #df$Mfactor <- as.factor(df$M)
  #df$intercept <- 1
  #jomomod <- jomo(Y = as.data.frame(df$Mfactor),
                  #X = as.data.frame(df[, c(10, 3, 4, 6)]),
                  #X2 = as.data.frame(df[, 3]),
                  #clus = as.data.frame(df[, 1]),
                  #nburn = 1000, nbetween = 1000, nimp = numimp, output = 0)
  #imp_mitml <- jomo2mitml.list(jomomod)
  #fit_i <- with(imp_mitml, geeglm(Y ~ A*df.Mfactor, family = "gaussian",
                                  #id = clus, corstr = "exchangeable"))
  #fit_jomo <- testEstimates(fit_i, extra.pars = T)
  
  #testmod <- jomoImpute(as.data.frame(df),
                        #formula = Mfactor ~ X*A*Y + (1 | cluster_ID))
  #imp_mitml <- mitmlComplete(testmod)
  #fit_i <- with(imp_mitml, geeglm(Y ~ A*Mfactor + X:A:Mfactor,
                                  #family = "gaussian",
                                  #id = cluster_ID, corstr = "exchangeable"))
  #fit_jomo <- testEstimates(fit_i, extra.pars = T)
  
  # Logit Gibbs samplers
  # Initialize and set priors
  numburn <- 1000
  thin <- 100
  numiter <- numburn + thin*numimp
  # Gamma hyperpriors
  c <- d <- 0.01
  
  # First do conditional on X + A + Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod4)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 4)
  
  # Initial values
  beta <- summary(impmod4)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y)
  ests12 <- array(NA, dim = c(numimp, 4))
  varests12 <- array(NA, dim = c(numimp, 4))

  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp12 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp12$M[is.na(dfimp12$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod12 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp12,
                      id = cluster_ID, corstr = "exchangeable")
      ests12[(h - numburn) / thin, ] <- coef(mod12)
      varests12[(h - numburn) / thin, ] <- (summary(mod12)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp12$M - 1/2) / omega
    
    # Update beta
    v <- solve(crossprod(sqrt(omega)*X) + T0)
    m <- v %*% (T0 %*% beta0 + t(sqrt(omega)*X) %*%
               (sqrt(omega)*(z - rep(b, size_clusters))))
    beta <- c(rmvnorm(1, m, v))
    
    # Update b
    vb <- 1 / (taub + tapply(omega, id, sum))
    mb <- vb*(tapply(omega*(z - X %*% beta), id, sum))
    b <- rnorm(num_clusters, mb, sqrt(vb))
    
    # Update taub
    taub <- rgamma(1, c + num_clusters/2, d + crossprod(b)/2)
    
  }
  
  # Next do conditional on X + A*Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod5)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 5)
  
  # Initial values
  beta <- summary(impmod5)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$A*df$Y)
  ests13 <- array(NA, dim = c(numimp, 4))
  varests13 <- array(NA, dim = c(numimp, 4))
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp13 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp13$M[is.na(dfimp13$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod13 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp13,
                      id = cluster_ID, corstr = "exchangeable")
      ests13[(h - numburn) / thin, ] <- coef(mod13)
      varests13[(h - numburn) / thin, ] <- (summary(mod13)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp13$M - 1/2) / omega
    
    # Update beta
    v <- solve(crossprod(sqrt(omega)*X) + T0)
    m <- v %*% (T0 %*% beta0 + t(sqrt(omega)*X) %*%
                  (sqrt(omega)*(z - rep(b, size_clusters))))
    beta <- c(rmvnorm(1, m, v))
    
    # Update b
    vb <- 1 / (taub + tapply(omega, id, sum))
    mb <- vb*(tapply(omega*(z - X %*% beta), id, sum))
    b <- rnorm(num_clusters, mb, sqrt(vb))
    
    # Update taub
    taub <- rgamma(1, c + num_clusters/2, d + crossprod(b)/2)
    
  }
  
  # Now do conditional on X*A + Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod6)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 5)
  
  # Initial values
  beta <- summary(impmod6)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$X*df$A)
  ests14 <- array(NA, dim = c(numimp, 4))
  varests14 <- array(NA, dim = c(numimp, 4))
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp14 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp14$M[is.na(dfimp14$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod14 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp14,
                      id = cluster_ID, corstr = "exchangeable")
      ests14[(h - numburn) / thin, ] <- coef(mod14)
      varests14[(h - numburn) / thin, ] <- (summary(mod14)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp14$M - 1/2) / omega
    
    # Update beta
    v <- solve(crossprod(sqrt(omega)*X) + T0)
    m <- v %*% (T0 %*% beta0 + t(sqrt(omega)*X) %*%
                  (sqrt(omega)*(z - rep(b, size_clusters))))
    beta <- c(rmvnorm(1, m, v))
    
    # Update b
    vb <- 1 / (taub + tapply(omega, id, sum))
    mb <- vb*(tapply(omega*(z - X %*% beta), id, sum))
    b <- rnorm(num_clusters, mb, sqrt(vb))
    
    # Update taub
    taub <- rgamma(1, c + num_clusters/2, d + crossprod(b)/2)
    
  }
  
  # Now do conditional on X*A + Y*A
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod7)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 6)
  
  # Initial values
  beta <- summary(impmod7)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$X*df$A, df$A*df$Y)
  ests15 <- array(NA, dim = c(numimp, 4))
  varests15 <- array(NA, dim = c(numimp, 4))
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp15 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp15$M[is.na(dfimp15$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod15 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp15,
                      id = cluster_ID, corstr = "exchangeable")
      ests15[(h - numburn) / thin, ] <- coef(mod15)
      varests15[(h - numburn) / thin, ] <- (summary(mod15)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp15$M - 1/2) / omega
    
    # Update beta
    v <- solve(crossprod(sqrt(omega)*X) + T0)
    m <- v %*% (T0 %*% beta0 + t(sqrt(omega)*X) %*%
                  (sqrt(omega)*(z - rep(b, size_clusters))))
    beta <- c(rmvnorm(1, m, v))
    
    # Update b
    vb <- 1 / (taub + tapply(omega, id, sum))
    mb <- vb*(tapply(omega*(z - X %*% beta), id, sum))
    b <- rnorm(num_clusters, mb, sqrt(vb))
    
    # Update taub
    taub <- rgamma(1, c + num_clusters/2, d + crossprod(b)/2)
    
  }
  
  # Now do conditional on X*A*Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod8)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 8)
  
  # Initial values
  beta <- summary(impmod8)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$X*df$A, df$X*df$Y,
             df$A*df$Y, df$X*df$A*df$Y)
  ests16 <- array(NA, dim = c(numimp, 4))
  varests16 <- array(NA, dim = c(numimp, 4))
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp16 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp16$M[is.na(dfimp16$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod16 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp16,
                      id = cluster_ID, corstr = "exchangeable")
      ests16[(h - numburn) / thin, ] <- coef(mod16)
      varests16[(h - numburn) / thin, ] <- (summary(mod16)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp16$M - 1/2) / omega
    
    # Update beta
    v <- solve(crossprod(sqrt(omega)*X) + T0)
    m <- v %*% (T0 %*% beta0 + t(sqrt(omega)*X) %*%
                  (sqrt(omega)*(z - rep(b, size_clusters))))
    beta <- c(rmvnorm(1, m, v))
    
    # Update b
    vb <- 1 / (taub + tapply(omega, id, sum))
    mb <- vb*(tapply(omega*(z - X %*% beta), id, sum))
    b <- rnorm(num_clusters, mb, sqrt(vb))
    
    # Update taub
    taub <- rgamma(1, c + num_clusters/2, d + crossprod(b)/2)
    
  }
  
  
  # Probit Gibbs sampler
  # Prior mean for beta
  #beta0 <- summary(impmod8)$coefficients[, 1]
  # Prior precision for beta
  #T0 <- diag(.01, 8)
  # Priors for Random Effects
  #nu0 <- 3				      # DF for IWish prior on Sigmab  
  #c0 <- 0.01			# Scale matrix for Wishart Prior on Sigmae
  
  # Initial values
  #beta <- summary(impmod8)$coefficients[, 1]
  #z <- rep(0, dim(df)[1])
  #b <- rep(0, num_clusters) # Random effects
  #taub <- 2 # Random effect precision
  
  # Save output
  #ests13 <- array(NA, dim = c(numimp, 4))
  #varests13 <- array(NA, dim = c(numimp, 4))
  
  # Posterior Var(beta) outside of loop
  #vbeta <- solve(T0 + crossprod(X, X))
  
  # Algorithm
  #for (h in 1:numiter) {
    
    # Impute missing data
    #dfimp13 <- df
    #est <- X %*% beta + rep(b, size_clusters)
    #estprob <- pnorm(est)
    #dfimp13$M[is.na(dfimp13$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    #if (h > numburn & h %% thin == 0) {
      #mod13 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp13,
                      #id = cluster_ID, corstr = "exchangeable")
      #ests13[(h - numburn) / thin, ] <- coef(mod13)
      #varests13[(h - numburn) / thin, ] <- (summary(mod13)$coefficients[, 2])^2
    #}
    
    # Draw latent variable z
    #muz <- X %*% beta + rep(b, size_clusters)
    #z[dfimp13$M == 0] <- qnorm(runif(dim(df)[1], 0, pnorm(0, muz)),
                               #muz)[dfimp13$M == 0]
    #z[dfimp13$M == 1] <- qnorm(runif(dim(df)[1], pnorm(0, muz), 1),
                               #muz)[dfimp13$M == 1]
    #z[z == -Inf] <- -100
    
    # Update beta
    #mbeta <- vbeta %*% (T0 %*% beta0 + crossprod(X, z - rep(b, size_clusters)))
    #beta <- c(rmvnorm(1, mbeta, vbeta))
    
    # Update b
    #taub12<-taub1/(1-rhob^2)  		                    # Prior precision of b1|b2
    #mb12<-rhob*sqrt(taub2/taub1)*b2 	                # Prior mean of b1|b2
    #vb1<-1/(nis+taub12)                              # Posterior var of b1|b2,y
    #mb1<-vb1*(taub12*mb12+tapply(z-X%*%beta,id,sum)) # Posterior mean of b1|b2,y
    #b1<-rnorm(n,mb1,sqrt(vb1))
    
    # Update b
    #vb <- 1 / (c(taub) + size_clusters)
    #vb <- 1 / (c(taub) + tapply(z, id, sum))
    #mb <- vb*(c(taub) + tapply((z - X %*% beta), id, sum))
    #mb <- vb*(tapply((z - X %*% beta), id, sum))
    #b <- rnorm(num_clusters, mb, sqrt(vb))
    
    # Update taub
    #Sigmab <- riwish(nu0 + num_clusters, c0 + crossprod(b))
    #taub <- 1 / Sigmab
    
  #}
  
  # Helper function to find t-value (>1.96 when num_clusters is small)
  find_tval <- function(ests, varests, numimp) {
    nu <- (numimp - 1)*(1 + numimp*colMeans(varests) /
                          (numimp + 1) / colVars(ests))^2
    nucom <- num_clusters - 2
    nuobs <- nucom*(nucom + 1) / (nucom + 3) /
      (1 + (numimp + 1)*colVars(ests) / numimp / colMeans(ests))
    nuadj <- ((1 / nu) + (1 / nuobs))^-1
    return(qt(0.975, nuadj))
  }
  
  # Output
  return(c(coef(mod1), coef(mod2), colMeans(ests3), colMeans(ests4),
           colMeans(ests5), colMeans(ests6), colMeans(ests7), colMeans(ests8),
           #mean(ests9),
           #colMeans(ests10[51:100, ]),
           #fit_jomo$estimates[1:4, 1],
           #fit_jomo$estimates[6, 1] - fit_jomo$estimates[5, 1],
           colMeans(ests12), colMeans(ests13), colMeans(ests14),
           colMeans(ests15), colMeans(ests16),
           summary(mod1)$coefficients[, 2],
           summary(mod2)$coefficients[, 2],
           sqrt(colVars(ests3) + (numimp+1)/(numimp-1)*colMeans(varests3)),
           sqrt(colVars(ests4) + (numimp+1)/(numimp-1)*colMeans(varests4)),
           sqrt(colVars(ests5) + (numimp+1)/(numimp-1)*colMeans(varests5)),
           sqrt(colVars(ests6) + (numimp+1)/(numimp-1)*colMeans(varests6)),
           sqrt(colVars(ests7) + (numimp+1)/(numimp-1)*colMeans(varests7)),
           sqrt(colVars(ests8) + (numimp+1)/(numimp-1)*colMeans(varests8)),
           #sqrt(var(ests9) + (numimp+1)/(numimp-1)*mean(varests9)),
           #sqrt(colVars(ests10[51:100, ]) + colMeans(varests10[51:100, ])),
           #fit_jomo$estimates[-5, 2],
           sqrt(colVars(ests12) + (numimp+1)/(numimp-1)*colMeans(varests12)),
           sqrt(colVars(ests13) + (numimp+1)/(numimp-1)*colMeans(varests13)),
           sqrt(colVars(ests14) + (numimp+1)/(numimp-1)*colMeans(varests14)),
           sqrt(colVars(ests15) + (numimp+1)/(numimp-1)*colMeans(varests15)),
           sqrt(colVars(ests16) + (numimp+1)/(numimp-1)*colMeans(varests16)),
           find_tval(ests3, varests3, numimp),
           find_tval(ests4, varests4, numimp),
           find_tval(ests5, varests5, numimp),
           find_tval(ests6, varests6, numimp),
           find_tval(ests7, varests7, numimp),
           find_tval(ests8, varests8, numimp),
           #find_tval(ests8, varests8, numimp),
           #find_tval(ests10[51:100, ], varests10[51:100, ], 50),
           #c(colQuantiles(ests12, probs = c(0.025, 0.975))),
           #find_tval(ests9, varests9),
           #quantile(ests10[21:200], c(.025, .975)),
           find_tval(ests12, varests12, numimp),
           find_tval(ests13, varests13, numimp),
           find_tval(ests14, varests13, numimp),
           find_tval(ests15, varests13, numimp),
           find_tval(ests16, varests13, numimp)
           ))
  
}


# Send simulations to computing cluster
nsims <- 500
ICC_out <- 0.1
ICC_mod <- 0.1
num_clusters <- 50
combos <- data.frame(trials = seq(1, nsims),
                     ICC_outs = rep(ICC_out, nsims),
                     ICC_mods = rep(ICC_mod, nsims),
                     num_clusterss = rep(num_clusters, nsims))
#i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
i <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
combo_i <- combos[(i), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, ICC_outs, ICC_mods,
                            num_clusterss))

# Output
#outfile <- paste("./Results/results_mod_XM_Iout_", ICC_out, "_Imod_",
                 #ICC_mod, "_nc_", num_clusters, "_",
                 #i, ".Rdata", sep = "")
outfile <-
  paste("/project/mharhaylab/blette/3_17_22/Results/results_out",
        "_Iout_", ICC_out, "_Imod_", ICC_mod, "_nc_", num_clusters,
        "_", i, ".Rdata", sep = "")
save(sim, file = outfile)