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
  Rlogit <- 1.2 + 0.5*df$X #- 0.2*df$Y
  Rprob <- exp(Rlogit) / (1 + exp(Rlogit))
  df$R <- rbinom(n, 1, Rprob)
  
  df$M <- NA
  df$M[df$R == 1] <- df$Mfull[df$R == 1]
  
  # Add interactions to data for JAV analyses
  df$AM <- df$A*df$M
  df$XAM <- df$X*df$A*df$M
  
  
  ###########################################################################
  # Data analysis
  # Truth
  truemod <- geeglm(Y ~ A*Mfull + X:A:Mfull, family = "gaussian", data = df,
                    id = cluster_ID, corstr = "exchangeable")
  
  ###################
  # Method 1: CRA-GEE
  mod1 <- geeglm(Y ~ A*M, family = "gaussian",
                 data = df[!is.na(df$M), ],
                 id = cluster_ID, corstr = "exchangeable")
  
  ################################
  # Methods 2-5: Single Imputation
  # Method 2: Single imputation under X + A + Y model
  impmod2 <- glm(M ~ X + A + Y, data = df[!is.na(df$M), ],
                 family = "binomial")
  dfimp2 <- df
  dfimp2$M[is.na(dfimp2$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod2, df[is.na(df$M), ],
                                              type = "response"))
  mod2 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp2,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 3: Single imputation under X + A*Y model
  impmod3 <- glm(M ~ X + A*Y, data = df[!is.na(df$M), ],
                 family = "binomial")
  dfimp3 <- df
  dfimp3$M[is.na(dfimp3$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod3, df[is.na(df$M), ],
                                              type = "response"))
  mod3 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp3,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 4: Single imputation under X*A + Y model
  impmod4 <- glm(M ~ X*A + Y, data = df[!is.na(df$M), ],
                 family = "binomial")
  dfimp4 <- df
  dfimp4$M[is.na(dfimp4$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod4, df[is.na(df$M), ],
                                              type = "response"))
  mod4 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp4,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 5: Single imputation under X*A + Y*A model
  impmod5 <- glm(M ~ X*A + Y*A, data = df[!is.na(df$M), ],
                 family = "binomial")
  dfimp5 <- df
  dfimp5$M[is.na(dfimp5$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod5, df[is.na(df$M), ],
                                              type = "response"))
  mod5 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp5,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 6: Single imputation under X*A*Y model
  impmod6 <- glm(M ~ X*A*Y, data = df[!is.na(df$M), ],
                 family = "binomial")
  dfimp6 <- df
  dfimp6$M[is.na(dfimp6$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod6, df[is.na(df$M), ],
                                              type = "response"))
  mod6 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp6,
                 id = cluster_ID, corstr = "exchangeable")
  
  
  # Method 3: Single Imputation - JAV
  #impmod3_AM <- glm(AM ~ X*A*Y, data = df[!is.na(df$M), ],
                    #family = "binomial")
  #impmod3_XAM <- glm(XAM ~ X*A*Y, data = df[!is.na(df$M), ],
                     #family = "gaussian")
  
  #dfimp3 <- df
  #dfimp3$M[is.na(dfimp3$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      #predict(impmod2_M, df[is.na(df$M), ],
                                              #type = "response"))
  #dfimp3$AM[is.na(dfimp3$AM)] <- rbinom(sum(is.na(df$AM)), 1,
                                      #predict(impmod3_AM, df[is.na(df$AM), ],
                                                #type = "response"))
  #dfimp3$XAM[is.na(dfimp3$XAM)] <- predict(impmod3_XAM, df[is.na(df$XAM), ])
  #mod3 <- geeglm(Y ~ A + M + AM, family = "gaussian", data = dfimp3,
                 #id = cluster_ID, corstr = "exchangeable")
  
  # Method 4: Single Imputation - JAV enhanced
  #impmod4_AM <- glm(AM ~ X*Y, data = df[!is.na(df$M) & df$A != 0, ],
                    #family = "binomial")
  
  #dfimp4 <- df
  #dfimp4$M[is.na(dfimp4$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      #predict(impmod2_M, df[is.na(df$M), ],
                                              #type = "response"))
  #dfimp4$AM[is.na(dfimp4$AM) & dfimp4$A == 0] <- 0
  #dfimp4$AM[is.na(dfimp4$AM) & dfimp4$A != 0] <-
    #rbinom(sum(is.na(df$AM) & df$A != 0), 1,
           #predict(impmod4_AM, df[is.na(df$AM) & df$A != 0, ],
                   #type = "response"))
  #mod4 <- geeglm(Y ~ A + M + AM, family = "gaussian", data = dfimp4,
                 #id = cluster_ID, corstr = "exchangeable")
  
  ###########################################
  # Methods 7-11: Multiple Imputation methods
  # Methods 12-16: Multi-level Multiple Imputation methods
  numimp <- 15
  ests7 <- ests8 <- ests9 <- ests10 <- ests11 <- ests12 <- ests13 <- ests14 <-
    ests15 <- ests16 <- ests17 <- ests18 <- ests19 <- ests20 <- ests21 <-
    array(NA, dim = c(numimp, 4))
  varests7 <- varests8 <- varests9 <- varests10 <- varests11 <- varests12 <-
    varests13 <- varests14 <- varests15 <- varests16 <- varests17 <-
    varests18 <- varests19 <- varests20 <- varests21 <-
    array(NA, dim = c(numimp, 4))
  
  # Define imputation models for the MMI methods
  # MI methods will use the same imputation models as for SI defined above
  impmod12 <- glmer(M ~ X + A + Y + (1 | cluster_ID), family = "binomial",
                    data = df[!is.na(df$M), ])
  impmod13 <- glmer(M ~ X + A*Y + (1 | cluster_ID), family = "binomial",
                    data = df[!is.na(df$M), ])
  impmod14 <- glmer(M ~ X*A + Y + (1 | cluster_ID), family = "binomial",
                    data = df[!is.na(df$M), ])
  impmod15 <- glmer(M ~ X*A + Y*A + (1 | cluster_ID), family = "binomial",
                    data = df[!is.na(df$M), ])
  impmod16 <- glmer(M ~ X*A*Y + (1 | cluster_ID), family = "binomial",
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
  #impmodbart <- bart(M ~ . - cluster_ID, df[!is.na(df$M), c(1, 3, 4, 6, 8)],
                     #keeptrees = TRUE)
  #predtest <- predict(impmodbart, df[is.na(df$M), c(1, 3, 4, 6, 8)],
                      #type = "ppd")
  
  # Run MI and MMI methods
  for (m in 1:numimp) {
    
    #############################
    # Multiple imputation methods
    # Method 7: Multiple Imputation under X + A + Y
    dfimp7 <- df
    dfimp7$M[is.na(dfimp7$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod2, df[is.na(df$M), ],
             type = "response"))
    mod7 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp7,
                   id = cluster_ID, corstr = "exchangeable")
    ests7[m, ] <- coef(mod7)
    varests7[m, ] <- (summary(mod7)$coefficients[, 2])^2
    
    # Method 8: Multiple Imputation under X + A*Y
    dfimp8 <- df
    dfimp8$M[is.na(dfimp8$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod3, df[is.na(df$M), ],
                                          type = "response"))
    mod8 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp8,
                   id = cluster_ID, corstr = "exchangeable")
    ests8[m, ] <- coef(mod8)
    varests8[m, ] <- (summary(mod8)$coefficients[, 2])^2
    
    # Method 9: Multiple Imputation under X*A + Y
    dfimp9 <- df
    dfimp9$M[is.na(dfimp9$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod4, df[is.na(df$M), ],
                                          type = "response"))
    mod9 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp9,
                   id = cluster_ID, corstr = "exchangeable")
    ests9[m, ] <- coef(mod9)
    varests9[m, ] <- (summary(mod9)$coefficients[, 2])^2
    
    # Method 10: Multiple Imputation under X*A + Y*A
    dfimp10 <- df
    dfimp10$M[is.na(dfimp10$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod5, df[is.na(df$M), ],
                                          type = "response"))
    mod10 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp10,
                    id = cluster_ID, corstr = "exchangeable")
    ests10[m, ] <- coef(mod10)
    varests10[m, ] <- (summary(mod10)$coefficients[, 2])^2
    
    # Method 11: Multiple Imputation under X*A*Y
    dfimp11 <- df
    dfimp11$M[is.na(dfimp11$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod6, df[is.na(df$M), ],
                                          type = "response"))
    mod11 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp11,
                    id = cluster_ID, corstr = "exchangeable")
    ests11[m, ] <- coef(mod11)
    varests11[m, ] <- (summary(mod11)$coefficients[, 2])^2
    
    ########################################
    # Multilevel multiple imputation methods
    # Method 12: Multilevel MI under X + A + Y
    dfimp12 <- df
    dfimp12$M[is.na(dfimp12$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod12, df[is.na(df$M), ],
                                          type = "response"))
    mod12 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp12,
                    id = cluster_ID, corstr = "exchangeable")
    ests12[m, ] <- coef(mod12)
    varests12[m, ] <- (summary(mod12)$coefficients[, 2])^2
    
    # Method 13: Multilevel MI under X + A*Y
    dfimp13 <- df
    dfimp13$M[is.na(dfimp13$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod13, df[is.na(df$M), ],
                                          type = "response"))
    mod13 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp13,
                    id = cluster_ID, corstr = "exchangeable")
    ests13[m, ] <- coef(mod13)
    varests13[m, ] <- (summary(mod13)$coefficients[, 2])^2
    
    # Method 14: Multilevel MI under X*A + Y
    dfimp14 <- df
    dfimp14$M[is.na(dfimp14$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod14, df[is.na(df$M), ],
                                          type = "response"))
    mod14 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp14,
                    id = cluster_ID, corstr = "exchangeable")
    ests14[m, ] <- coef(mod14)
    varests14[m, ] <- (summary(mod14)$coefficients[, 2])^2
    
    # Method 15: Multilevel MI under X*A + Y*A
    dfimp15 <- df
    dfimp15$M[is.na(dfimp15$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod15, df[is.na(df$M), ],
                                          type = "response"))
    mod15 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp15,
                    id = cluster_ID, corstr = "exchangeable")
    ests15[m, ] <- coef(mod15)
    varests15[m, ] <- (summary(mod15)$coefficients[, 2])^2
    
    # Method 16: Multilevel MI under X*A*Y
    dfimp16 <- df
    dfimp16$M[is.na(dfimp16$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod16, df[is.na(df$M), ],
                                          type = "response"))
    mod16 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp16,
                    id = cluster_ID, corstr = "exchangeable")
    ests16[m, ] <- coef(mod16)
    varests16[m, ] <- (summary(mod16)$coefficients[, 2])^2
    
    # Method 9: BART imputation
    #dfimp9 <- df
    #dfimp9$M[is.na(dfimp9$M)] <-
      #rbinom(sum(is.na(df$M)), 1, colMeans(predtest))
    #mod9 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp9,
                   #id = cluster_ID, corstr = "exchangeable")
    #ests9[m, ] <- coef(mod9)
    #varests9[m, ] <- (summary(mod9)$coefficients[, 2])^2
    
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
  
  
  #################################################
  # Bayesian MMI methods using logit Gibbs samplers
  # Initialize and set priors
  numburn <- 1000
  thin <- 100
  numiter <- numburn + thin*numimp
  # Gamma hyperpriors
  c <- d <- 0.01
  
  # Method 17: BMMI under X + A + Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod12)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 4)
  
  # Initial values
  beta <- summary(impmod12)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y)

  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp17 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp17$M[is.na(dfimp17$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod17 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp17,
                      id = cluster_ID, corstr = "exchangeable")
      ests17[(h - numburn) / thin, ] <- coef(mod17)
      varests17[(h - numburn) / thin, ] <- (summary(mod17)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp17$M - 1/2) / omega
    
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
  
  # Method 18: BMMI under X + A*Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod13)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 5)
  
  # Initial values
  beta <- summary(impmod13)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$A*df$Y)
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp18 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp18$M[is.na(dfimp18$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod18 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp18,
                      id = cluster_ID, corstr = "exchangeable")
      ests18[(h - numburn) / thin, ] <- coef(mod18)
      varests18[(h - numburn) / thin, ] <- (summary(mod18)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp18$M - 1/2) / omega
    
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
  
  # Method 19: BMMI under X*A + Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod14)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 5)
  
  # Initial values
  beta <- summary(impmod14)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$X*df$A)
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp19 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp19$M[is.na(dfimp19$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod19 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp19,
                      id = cluster_ID, corstr = "exchangeable")
      ests19[(h - numburn) / thin, ] <- coef(mod19)
      varests19[(h - numburn) / thin, ] <- (summary(mod19)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp19$M - 1/2) / omega
    
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
  
  # Method 20: BMMI under X*A + Y*A
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod15)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 6)
  
  # Initial values
  beta <- summary(impmod15)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$X*df$A, df$A*df$Y)
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp20 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp20$M[is.na(dfimp20$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod20 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp20,
                      id = cluster_ID, corstr = "exchangeable")
      ests20[(h - numburn) / thin, ] <- coef(mod20)
      varests20[(h - numburn) / thin, ] <- (summary(mod20)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp20$M - 1/2) / omega
    
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
  
  # Method 21: BMMI under X*A*Y
  # Set priors
  # Prior mean for beta
  beta0 <- summary(impmod16)$coefficients[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 8)
  
  # Initial values
  beta <- summary(impmod16)$coefficients[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- df$cluster_ID
  X <- cbind(1, df$X, df$A, df$Y, df$X*df$A, df$X*df$Y,
             df$A*df$Y, df$X*df$A*df$Y)
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    dfimp21 <- df
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(df$M)]) / (1 + exp(est[is.na(df$M)]))
    dfimp21$M[is.na(dfimp21$M)] <- rbinom(sum(is.na(df$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod21 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp21,
                      id = cluster_ID, corstr = "exchangeable")
      ests21[(h - numburn) / thin, ] <- coef(mod21)
      varests21[(h - numburn) / thin, ] <- (summary(mod21)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(df)[1], 1, mu)
    z <- (dfimp21$M - 1/2) / omega
    
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
  return(c(coef(mod1), coef(mod2), coef(mod3), coef(mod4), coef(mod5),
           coef(mod6), colMeans(ests7), colMeans(ests8), colMeans(ests9),
           colMeans(ests10), colMeans(ests11), colMeans(ests12),
           colMeans(ests13), colMeans(ests14), colMeans(ests15),
           colMeans(ests16), colMeans(ests17), colMeans(ests18),
           colMeans(ests19), colMeans(ests20), colMeans(ests21),
           summary(mod1)$coefficients[, 2], summary(mod2)$coefficients[, 2],
           summary(mod3)$coefficients[, 2], summary(mod4)$coefficients[, 2],
           summary(mod5)$coefficients[, 2], summary(mod6)$coefficients[, 2],
           sqrt(colVars(ests7) + (numimp+1)/(numimp-1)*colMeans(varests7)),
           sqrt(colVars(ests8) + (numimp+1)/(numimp-1)*colMeans(varests8)),
           sqrt(colVars(ests9) + (numimp+1)/(numimp-1)*colMeans(varests9)),
           sqrt(colVars(ests10) + (numimp+1)/(numimp-1)*colMeans(varests10)),
           sqrt(colVars(ests11) + (numimp+1)/(numimp-1)*colMeans(varests11)),
           sqrt(colVars(ests12) + (numimp+1)/(numimp-1)*colMeans(varests12)),
           sqrt(colVars(ests13) + (numimp+1)/(numimp-1)*colMeans(varests13)),
           sqrt(colVars(ests14) + (numimp+1)/(numimp-1)*colMeans(varests14)),
           sqrt(colVars(ests15) + (numimp+1)/(numimp-1)*colMeans(varests15)),
           sqrt(colVars(ests16) + (numimp+1)/(numimp-1)*colMeans(varests16)),
           sqrt(colVars(ests17) + (numimp+1)/(numimp-1)*colMeans(varests17)),
           sqrt(colVars(ests18) + (numimp+1)/(numimp-1)*colMeans(varests18)),
           sqrt(colVars(ests19) + (numimp+1)/(numimp-1)*colMeans(varests19)),
           sqrt(colVars(ests20) + (numimp+1)/(numimp-1)*colMeans(varests20)),
           sqrt(colVars(ests21) + (numimp+1)/(numimp-1)*colMeans(varests21)),
           find_tval(ests7, varests7, numimp),
           find_tval(ests8, varests8, numimp),
           find_tval(ests9, varests9, numimp),
           find_tval(ests10, varests10, numimp),
           find_tval(ests11, varests11, numimp),
           find_tval(ests12, varests12, numimp),
           find_tval(ests13, varests13, numimp),
           find_tval(ests14, varests14, numimp),
           find_tval(ests15, varests15, numimp),
           find_tval(ests16, varests16, numimp),
           find_tval(ests17, varests17, numimp),
           find_tval(ests18, varests18, numimp),
           find_tval(ests19, varests19, numimp),
           find_tval(ests20, varests20, numimp),
           find_tval(ests21, varests21, numimp)
           ))
  
}


# Send simulations to computing cluster
nsims <- 1000
ICC_out <- 0.1
ICC_mod <- 0.1
num_clusters <- 100
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
  paste("/project/mharhaylab/blette/4_7_22/Results/results_notout",
        "_Iout_", ICC_out, "_Imod_", ICC_mod, "_nc_", num_clusters,
        "_", i, ".Rdata", sep = "")
save(sim, file = outfile)