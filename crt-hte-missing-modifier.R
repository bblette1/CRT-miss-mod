# Exploratory simulation of HTE in CRT w/ missing moderator
# Load libraries
library(dplyr)
library(geepack)
library(ggplot2)
library(lme4)
#library(mice)

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
    0.7*df$X*df$Mfull*df$A +
    rep(alpha1, size_clusters) + rnorm(n, 0, sqrt(Y_resvar))
  
  # Missingness
  Rlogit <- 0.5 + 0.5*df$X
  Rprob <- exp(Rlogit) / (1 + exp(Rlogit))
  df$R <- rbinom(n, 1, Rprob)
  
  df$M <- NA
  df$M[df$R == 1] <- df$Mfull[df$R == 1]
  
  
  # Data analysis
  # Truth
  truemod <- geeglm(Y ~ A*Mfull, family = "gaussian", data = df,
                    id = cluster_ID, corstr = "exchangeable")
  
  # Method 1: CRA-GEE
  mod1 <- geeglm(Y ~ A*M, family = "gaussian",
                 data = df[!is.na(df$M), ],
                 id = cluster_ID, corstr = "exchangeable")
  
  # Method 2: Single Imputation
  impmod2 <- glm(M ~ X*A*Y, data = df[!is.na(df$M), ],
                 family = "binomial")
  dfimp2 <- df
  dfimp2$M[is.na(dfimp2$M)] <- rbinom(sum(is.na(df$M)), 1,
                                      predict(impmod2, df[is.na(df$M), ],
                                              type = "response"))
  mod2 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp2,
                 id = cluster_ID, corstr = "exchangeable")
  
  # Methods 3-8: Multiple Imputation methods
  numimp <- 15
  ests3 <- ests4 <- ests5 <- ests6 <- ests7 <- ests8 <- rep(NA, numimp)
  varests3 <- varests4 <- varests5 <- varests6 <-
    varests7 <- varests8 <-rep(NA, numimp)
  
  impmod4 <- glmer(M ~ X + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  #impmod3 <- geeglm(M ~ X + A*Y, family = binomial(),
  #data = df[!is.na(df$M), ],
  #id = cluster_ID, corstr = "exchangeable")
  impmod5 <- glmer(M ~ X + Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  impmod6 <- glmer(M ~ X + A + Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  impmod7 <- glmer(M ~ X + A*Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  impmod8 <- glmer(M ~ X*A*Y + (1 | cluster_ID), family = "binomial",
                   data = df[!is.na(df$M), ])
  
  for (m in 1:numimp) {
    
    # Method 3: Multiple Imputation
    dfimp3 <- df
    dfimp3$M[is.na(dfimp3$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod2, df[is.na(df$M), ],
             type = "response"))
    mod3 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp3,
                   id = cluster_ID, corstr = "exchangeable")
    ests3[m] <- coef(mod3)[4]
    varests3[m] <- (summary(mod3)$coefficients[4, 2])^2
    
    # Method 4: Multilevel MI cond. on X
    dfimp4 <- df
    dfimp4$M[is.na(dfimp4$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod4, df[is.na(df$M), ],
                                          type = "response"))
    mod4 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp4,
                   id = cluster_ID, corstr = "exchangeable")
    ests4[m] <- coef(mod4)[4]
    varests4[m] <- (summary(mod4)$coefficients[4, 2])^2
    
    # Method 5: Multilevel MI cond. on X and Y
    dfimp5 <- df
    dfimp5$M[is.na(dfimp5$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod5, df[is.na(df$M), ],
                                          type = "response"))
    mod5 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp5,
                   id = cluster_ID, corstr = "exchangeable")
    ests5[m] <- coef(mod5)[4]
    varests5[m] <- (summary(mod5)$coefficients[4, 2])^2
    
    # Method 6: Multilevel MI cond. on X, Y, and A
    dfimp6 <- df
    dfimp6$M[is.na(dfimp6$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod6, df[is.na(df$M), ],
                                          type = "response"))
    mod6 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp6,
                   id = cluster_ID, corstr = "exchangeable")
    ests6[m] <- coef(mod6)[4]
    varests6[m] <- (summary(mod6)$coefficients[4, 2])^2
    
    # Method 7: Multilevel MI cond. on X and A*Y
    dfimp7 <- df
    dfimp7$M[is.na(dfimp7$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod7, df[is.na(df$M), ],
                                          type = "response"))
    mod7 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp7,
                   id = cluster_ID, corstr = "exchangeable")
    ests7[m] <- coef(mod7)[4]
    varests7[m] <- (summary(mod7)$coefficients[4, 2])^2
    
    # Method 8: Multilevel MI cond. on X*Y*A
    dfimp8 <- df
    dfimp8$M[is.na(dfimp8$M)] <-
      rbinom(sum(is.na(df$M)), 1, predict(impmod8, df[is.na(df$M), ],
                                          type = "response"))
    mod8 <- geeglm(Y ~ A*M, family = "gaussian", data = dfimp8,
                   id = cluster_ID, corstr = "exchangeable")
    ests8[m] <- coef(mod8)[4]
    varests8[m] <- (summary(mod8)$coefficients[4, 2])^2
    
  }
  
  find_tval <- function(ests, varests) {
    nu <- (numimp - 1)*(1 + numimp*mean(varests) /
                          (numimp + 1) / var(ests))^2
    nucom <- num_clusters - 2
    nuobs <- nucom*(nucom + 1) / (nucom + 3) /
      (1 + (numimp + 1)*var(ests) / numimp / mean(ests))
    nuadj <- ((1 / nu) + (1 / nuobs))^-1
    return(qt(0.975, nuadj))
  }
  
  # Output
  return(c(coef(mod1)[4], coef(mod2)[4], mean(ests3), mean(ests4),
           mean(ests5), mean(ests6), mean(ests7), mean(ests8),
           summary(mod1)$coefficients[4, 2],
           summary(mod2)$coefficients[4, 2],
           sqrt(var(ests3) + (numimp+1)/(numimp-1)*mean(varests3)),
           sqrt(var(ests4) + (numimp+1)/(numimp-1)*mean(varests4)),
           sqrt(var(ests5) + (numimp+1)/(numimp-1)*mean(varests5)),
           sqrt(var(ests6) + (numimp+1)/(numimp-1)*mean(varests6)),
           sqrt(var(ests7) + (numimp+1)/(numimp-1)*mean(varests7)),
           sqrt(var(ests8) + (numimp+1)/(numimp-1)*mean(varests8)),
           find_tval(ests3, varests3), find_tval(ests4, varests4),
           find_tval(ests5, varests5), find_tval(ests6, varests6),
           find_tval(ests7, varests7), find_tval(ests8, varests8)))
  
}


# Send simulations to computing cluster
nsims <- 10
ICC_out <- 0.1
ICC_mod <- 0.1
num_clusters <- 20
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
#outfile <- paste("./Results/results_mod_Iout_", ICC_out, "_Imod_",
                 #ICC_mod, "_nc_", num_clusters, "_",
                 #i, ".Rdata", sep = "")
outfile <-
  paste("/project/mharhaylab/blette/1_20_22/Results/results_mod_Iout_",
        ICC_out, "_Imod_", ICC_mod, "_nc_", num_clusters, "_",
        i, ".Rdata", sep = "")
save(sim, file = outfile)