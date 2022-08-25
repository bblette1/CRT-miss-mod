# Load libraries
rm(list = ls())
library(dplyr)
library(geepack)
library(GLMMadaptive)
library(lme4)
library(matrixStats)

# Helper function to find t-value for MI (>1.96 when num_clusters is small)
find_tval <- function(ests, varests, numimp, num_clusters) {
  
  nu <- (numimp - 1)*(1 + numimp*colMeans(varests) /
                        (numimp + 1) / colVars(ests))^2
  nucom <- num_clusters - 2
  nuobs <- nucom*(nucom + 1) / (nucom + 3) /
    (1 + (numimp + 1)*colVars(ests) / numimp / colMeans(ests))
  nuadj <- ((1 / nu) + (1 / nuobs))^-1
  
  outvec <- suppressWarnings(qt(0.975, nuadj))
  outvec[is.na(outvec)] <- mean(suppressWarnings(qt(0.975, nuadj)), na.rm = T)
  return(outvec)

}

# Load data
load("/Users/blette/Downloads/ICPSR_36158/DS0002/36158-0002-Data.rda")

# Set seed
set.seed(82322)

# Data wrangling
levels(da36158.0002$CONDITION) <- levels(da36158.0002$CONDITION)[2:1]
dat_long <- da36158.0002[da36158.0002$WAVE %in% c("(1) BASELINE",
                                                  "(3) 12-MONTH FOLLOW-UP"), ]
dat_long$wavenum <- NA
dat_long$wavenum[dat_long$WAVE == "(1) BASELINE"] <- 1
dat_long$wavenum[dat_long$WAVE == "(3) 12-MONTH FOLLOW-UP"] <- 2
dat_wide <- reshape(dat_long, idvar = "ADMINLINK", timevar = "wavenum",
                    direction = "wide")
dat_wide$sleepdiff <- dat_wide$CVPH_BEDHRS.2 - dat_wide$CVPH_BEDHRS.1
dat_wide$timediff <- dat_wide$SCWM_TIMEALLI.2 - dat_wide$SCWM_TIMEALLI.1

mod <- geeglm(timediff ~ CONDITION.1*MANAGER.1, family = "gaussian",
              data = dat_wide, id = RMZCLUSTER.1, corstr = "exchangeable")

# Make small data set for illustration
# WM_JSTR1R.1 is decision authority and WM_CWH2R.1 is job autonomy
dat <- dat_wide %>%
  dplyr::select(RMZCLUSTER.1, CONDITION.1, MANAGER.1, timediff, WM_JSTR1R.1,
                WM_CWH2R.1)
colnames(dat) <- c("Cluster", "A", "M", "Y", "X1", "X2")
dat <- dat[complete.cases(dat), ]
dat$M <- as.numeric(substr(dat$M, 2, 2))
dat$A_bin <- 0
dat$A_bin[dat$A == "(1) INTERVENTION"] <- 1
dat <- dat[order(dat$Cluster, decreasing = FALSE), ]

# Model with full data
mod <- geeglm(Y ~ A*M, family = "gaussian",
              data = dat, id = Cluster, corstr = "exchangeable")


#############################################################################
# Scenario 1: MCAR
nsims <- 100
cca_ests1 <- rep(NA, nsims)
cca_lower1 <- rep(NA, nsims)
cca_upper1 <- rep(NA, nsims)
si_ests1 <- rep(NA, nsims)
si_lower1 <- rep(NA, nsims)
si_upper1 <- rep(NA, nsims)
mi_ests1 <- rep(NA, nsims)
mi_lower1 <- rep(NA, nsims)
mi_upper1 <- rep(NA, nsims)
mmi_ests1 <- rep(NA, nsims)
mmi_lower1 <- rep(NA, nsims)
mmi_upper1 <- rep(NA, nsims)
bmmi_ests1 <- rep(NA, nsims)
bmmi_lower1 <- rep(NA, nsims)
bmmi_upper1 <- rep(NA, nsims)

for (i in 1:nsims) {
  
  dat1 <- dat
  misslist <- sample(1:dim(dat1)[1], dim(dat1)[1]/5, replace = FALSE)
  dat1$M[misslist] <- NA
  
  # CCA
  mod_cca <- geeglm(Y ~ A*M, family = "gaussian",
                    data = dat1[complete.cases(dat1), ],
                    id = Cluster, corstr = "exchangeable")
  cca_ests1[i] <- mod_cca$coefficients[4]
  cca_lower1[i] <- mod_cca$coefficients[4] -
                   1.96*summary(mod_cca)$coefficients[4, 2]
  cca_upper1[i] <- mod_cca$coefficients[4] +
                   1.96*summary(mod_cca)$coefficients[4, 2]
  
  # Single imputation
  impmod <- glm(M ~ Y*A, family = "binomial", data = dat1[!is.na(dat1$M), ])
  datimp1 <- dat1
  datimp1$M[is.na(datimp1$M)] <- rbinom(length(datimp1$M[is.na(datimp1$M)]), 1,
                                        predict(impmod, dat1[is.na(dat1$M), ],
                                                type = "response"))
  mod_si <- geeglm(Y ~ A*M, family = "gaussian",
                   data = datimp1,
                   id = Cluster, corstr = "exchangeable")
  si_ests1[i] <- mod_si$coefficients[4]
  si_lower1[i] <- mod_si$coefficients[4] -
                  1.96*summary(mod_si)$coefficients[4, 2]
  si_upper1[i] <- mod_si$coefficients[4] +
                  1.96*summary(mod_si)$coefficients[4, 2]
  
  # Multiple imputation
  numimp <- 15
  ests2 <- array(NA, dim = c(numimp, 4))
  varests2 <- array(NA, dim = c(numimp, 4))
  
  for (m in 1:numimp) {
    
    datimp2 <- dat1
    datimp2$M[is.na(datimp2$M)] <-
      rbinom(length(datimp2$M[is.na(datimp2$M)]), 1,
                    predict(impmod, dat1[is.na(dat1$M), ],
                                                  type = "response"))
    mod_mi <- geeglm(Y ~ A*M, family = "gaussian",
                     data = datimp2,
                     id = Cluster, corstr = "exchangeable")
    ests2[m, ] <- coef(mod_mi)
    varests2[m, ] <- (summary(mod_mi)$coefficients[, 2])^2
    
  }
  
  mi_ests1[i] <- colMeans(ests2)[4]
  mi_lowervec <- colMeans(ests2) -
    sqrt(colVars(ests2) + (numimp+1)/(numimp-1)*colMeans(varests2))*
    find_tval(ests2, varests2, numimp, length(unique(dat$Cluster)))
  mi_lower1[i] <- mi_lowervec[4]
  mi_uppervec <- colMeans(ests2) +
    sqrt(colVars(ests2) + (numimp+1)/(numimp-1)*colMeans(varests2))*
    find_tval(ests2, varests2, numimp, length(unique(dat$Cluster)))
  mi_upper1[i] <- mi_uppervec[4]
  
  # Multilevel multiple imputation
  numimp <- 15
  ests3 <- array(NA, dim = c(numimp, 4))
  varests3 <- array(NA, dim = c(numimp, 4))
  
  mmi_impmod <- mixed_model(M ~ Y*A_bin, random = ~ 1 | Cluster,
                            family = binomial(), 
                            data = dat1[!is.na(dat1$M), ], iter_EM = 100)
  
  for (m in 1:numimp) {
    
    datimp3 <- dat1
    datimp3$M[is.na(datimp3$M)] <-
      rbinom(length(datimp3$M[is.na(datimp3$M)]), 1,
                    predict(mmi_impmod, dat1[is.na(dat1$M), ]))
    mod_mmi <- geeglm(Y ~ A*M, family = "gaussian",
                      data = datimp3,
                      id = Cluster, corstr = "exchangeable")
    ests3[m, ] <- coef(mod_mmi)
    varests3[m, ] <- (summary(mod_mmi)$coefficients[, 2])^2
    
  }
  
  mmi_ests1[i] <- colMeans(ests3)[4]
  mmi_lowervec <- colMeans(ests3) -
    sqrt(colVars(ests3) + (numimp+1)/(numimp-1)*colMeans(varests3))*
    find_tval(ests3, varests3, numimp, length(unique(dat$Cluster)))
  mmi_lower1[i] <- mmi_lowervec[4]
  mmi_uppervec <- colMeans(ests3) +
    sqrt(colVars(ests3) + (numimp+1)/(numimp-1)*colMeans(varests3))*
    find_tval(ests3, varests3, numimp, length(unique(dat$Cluster)))
  mmi_upper1[i] <- mmi_uppervec[4]
  
  # Bayesian multilevel multiple imputation
  ests4 <- array(NA, dim = c(numimp, 4))
  varests4 <- array(NA, dim = c(numimp, 4))
  num_clusters <- length(unique(dat$Cluster))
  size_clusters <- table(dat$Cluster)
  
  # Initialize and set priors
  numburn <- 1000
  thin <- 100
  numiter <- numburn + thin*numimp
  # Gamma hyperpriors
  c <- d <- 0.01
  
  # Set priors
  # Prior mean for beta
  beta0 <- summary(mmi_impmod)$coef_table[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 4)
  
  # Initial values
  beta <- summary(mmi_impmod)$coef_table[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- dat$Cluster
  X <- cbind(1, dat$Y, dat$A_bin, dat$A_bin*dat$Y)
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    datimp4 <- dat1
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(dat1$M)]) / (1 + exp(est[is.na(dat1$M)]))
    datimp4$M[is.na(datimp4$M)] <- rbinom(sum(is.na(dat1$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod_bmmi <- geeglm(Y ~ A_bin*M, family = "gaussian", data = datimp4,
                         id = Cluster, corstr = "exchangeable")
      ests4[(h - numburn) / thin, ] <- coef(mod_bmmi)
      varests4[(h - numburn) / thin, ] <-
        (summary(mod_bmmi)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(dat1)[1], 1, mu)
    z <- (datimp4$M - 1/2) / omega
    
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
  
  bmmi_ests1[i] <- colMeans(ests4)[4]
  bmmi_lowervec <- colMeans(ests4) -
    sqrt(colVars(ests4) + (numimp+1)/(numimp-1)*colMeans(varests4))*
    find_tval(ests4, varests4, numimp, length(unique(dat$Cluster)))
  bmmi_lower1[i] <- bmmi_lowervec[4]
  bmmi_uppervec <- colMeans(ests4) +
    sqrt(colVars(ests4) + (numimp+1)/(numimp-1)*colMeans(varests4))*
    find_tval(ests4, varests4, numimp, length(unique(dat$Cluster)))
  bmmi_upper1[i] <- bmmi_uppervec[4]
  
}

plotdat <- data.frame(type = rep("Average"),
                      scenario = rep("Scenario 1", 6),
                      method = c("Real Data", "CCA", "SI", "MI", "MMI",
                                 "B-MMI"),
                      pointest = c(mod$coefficients[4], mean(cca_ests1),
                                   mean(si_ests1), mean(mi_ests1),
                                   mean(mmi_ests1), mean(bmmi_ests1)),
                      lower = c(mod$coefficients[4] -
                                  1.96*summary(mod)$coefficients[4, 2],
                                mean(cca_lower1),
                                mean(si_lower1), mean(mi_lower1),
                                mean(mmi_lower1),
                                mean(bmmi_lower1)),
                      upper = c(mod$coefficients[4] +
                                  1.96*summary(mod)$coefficients[4, 2],
                                mean(cca_upper1),
                                mean(si_upper1), mean(mi_upper1),
                                mean(mmi_upper1),
                                mean(bmmi_upper1)))

plotdat_single <- data.frame(type = rep("Single Draw"),
                      scenario = rep("Scenario 1", 6),
                      method = c("Real Data", "CCA", "SI", "MI", "MMI",
                                 "B-MMI"),
                      pointest = c(mod$coefficients[4], cca_ests1[1],
                                   si_ests1[1], mi_ests1[1],
                                   mmi_ests1[1], bmmi_ests1[1]),
                      lower = c(mod$coefficients[4] -
                                  1.96*summary(mod)$coefficients[4, 2],
                                cca_lower1[1],
                                si_lower1[1], mi_lower1[1],
                                mmi_lower1[1],
                                bmmi_lower1[1]),
                      upper = c(mod$coefficients[4] +
                                  1.96*summary(mod)$coefficients[4, 2],
                                cca_upper1[1],
                                si_upper1[1], mi_upper1[1],
                                mmi_upper1[1],
                                bmmi_upper1[1]))

ggplot(data = plotdat, aes(y = method, x = pointest, xmin = lower,
                          xmax = upper)) +
  geom_point() + 
  geom_errorbarh(height = .1) +
  scale_y_discrete(limits = rev(plotdat$method)) +
  labs(x = 'Interaction coefficient', y = 'Method') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5)

ggplot(data = plotdat_single, aes(y = method, x = pointest, xmin = lower,
                                  xmax = upper)) +
  geom_point() + 
  geom_errorbarh(height = .1) +
  scale_y_discrete(limits = rev(plotdat$method)) +
  labs(x = 'Interaction coefficient', y = 'Method') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5)


#############################################################################
# Scenario 2: Simple MAR
dat$X1_bin <- 0
dat$X1_bin[dat$X1 == "(5) STRONGLY AGREE" | dat$X1 == "(4) AGREE"] <- 1
dat$X2_bin <- 0
dat$X2_bin[dat$X2 == "(5) VERY MUCH" | dat$X2 == "(4) MUCH"] <- 1

cca_ests2 <- rep(NA, nsims)
cca_lower2 <- rep(NA, nsims)
cca_upper2 <- rep(NA, nsims)
si_ests2 <- rep(NA, nsims)
si_lower2 <- rep(NA, nsims)
si_upper2 <- rep(NA, nsims)
mi_ests2 <- rep(NA, nsims)
mi_lower2 <- rep(NA, nsims)
mi_upper2 <- rep(NA, nsims)
mmi_ests2 <- rep(NA, nsims)
mmi_lower2 <- rep(NA, nsims)
mmi_upper2 <- rep(NA, nsims)
bmmi_ests2 <- rep(NA, nsims)
bmmi_lower2 <- rep(NA, nsims)
bmmi_upper2 <- rep(NA, nsims)

for (i in 1:nsims) {
  
  dat2 <- dat
  Rlogit <- 2 + 0.5*dat2$Y - 0.6*dat2$X1_bin - 0.3*dat2$X2_bin
  Rprob <- exp(Rlogit) / (1 + exp(Rlogit))
  dat2$R <- rbinom(dim(dat2)[1], 1, Rprob)
  dat2$M[dat2$R == 0] <- NA
  
  # CCA
  mod_cca <- geeglm(Y ~ A*M, family = "gaussian",
                    data = dat2[complete.cases(dat2), ],
                    id = Cluster, corstr = "exchangeable")
  cca_ests2[i] <- mod_cca$coefficients[4]
  cca_lower2[i] <- mod_cca$coefficients[4] -
    1.96*summary(mod_cca)$coefficients[4, 2]
  cca_upper2[i] <- mod_cca$coefficients[4] +
    1.96*summary(mod_cca)$coefficients[4, 2]
  
  # Single imputation
  impmod <- glm(M ~ Y*A*X1_bin + Y*A*X2_bin, family = "binomial",
                data = dat2[!is.na(dat2$M), ])
  datimp5 <- dat2
  datimp5$M[is.na(datimp5$M)] <- rbinom(length(datimp5$M[is.na(datimp5$M)]), 1,
                                        predict(impmod, dat2[is.na(dat2$M), ],
                                                type = "response"))
  mod_si <- geeglm(Y ~ A*M, family = "gaussian",
                   data = datimp5,
                   id = Cluster, corstr = "exchangeable")
  si_ests2[i] <- mod_si$coefficients[4]
  si_lower2[i] <- mod_si$coefficients[4] -
    1.96*summary(mod_si)$coefficients[4, 2]
  si_upper2[i] <- mod_si$coefficients[4] +
    1.96*summary(mod_si)$coefficients[4, 2]
  
  # Multiple imputation
  numimp <- 15
  ests6 <- array(NA, dim = c(numimp, 4))
  varests6 <- array(NA, dim = c(numimp, 4))
  
  for (m in 1:numimp) {
    
    datimp6 <- dat2
    datimp6$M[is.na(datimp6$M)] <-
      rbinom(length(datimp6$M[is.na(datimp6$M)]), 1,
                              predict(impmod, dat2[is.na(dat2$M), ],
                                                  type = "response"))
    mod_mi <- geeglm(Y ~ A*M, family = "gaussian",
                     data = datimp6,
                     id = Cluster, corstr = "exchangeable")
    ests6[m, ] <- coef(mod_mi)
    varests6[m, ] <- (summary(mod_mi)$coefficients[, 2])^2
    
  }
  
  mi_ests2[i] <- colMeans(ests6)[4]
  mi_lowervec <- colMeans(ests6) -
    sqrt(colVars(ests6) + (numimp+1)/(numimp-1)*colMeans(varests6))*
    find_tval(ests6, varests6, numimp, length(unique(dat$Cluster)))
  mi_lower2[i] <- mi_lowervec[4]
  mi_uppervec <- colMeans(ests6) +
    sqrt(colVars(ests6) + (numimp+1)/(numimp-1)*colMeans(varests6))*
    find_tval(ests6, varests6, numimp, length(unique(dat$Cluster)))
  mi_upper2[i] <- mi_uppervec[4]
  
  # Multilevel multiple imputation
  numimp <- 15
  ests7 <- array(NA, dim = c(numimp, 4))
  varests7 <- array(NA, dim = c(numimp, 4))
  
  mmi_impmod <- mixed_model(M ~ Y*A*X1_bin + Y*A*X2_bin, random = ~ 1 | Cluster,
                            family = binomial(), 
                            data = dat2[!is.na(dat2$M), ], iter_EM = 100)
  
  for (m in 1:numimp) {
    
    datimp7 <- dat2
    datimp7$M[is.na(datimp7$M)] <-
      rbinom(length(datimp7$M[is.na(datimp7$M)]), 1,
                              predict(mmi_impmod, dat2[is.na(dat2$M), ]))
    mod_mmi <- geeglm(Y ~ A*M, family = "gaussian",
                      data = datimp7,
                      id = Cluster, corstr = "exchangeable")
    ests7[m, ] <- coef(mod_mmi)
    varests7[m, ] <- (summary(mod_mmi)$coefficients[, 2])^2
    
  }
  
  mmi_ests2[i] <- colMeans(ests7)[4]
  mmi_lowervec <- colMeans(ests7) -
    sqrt(colVars(ests7) + (numimp+1)/(numimp-1)*colMeans(varests7))*
    find_tval(ests7, varests7, numimp, length(unique(dat$Cluster)))
  mmi_lower2[i] <- mmi_lowervec[4]
  mmi_uppervec <- colMeans(ests7) +
    sqrt(colVars(ests7) + (numimp+1)/(numimp-1)*colMeans(varests7))*
    find_tval(ests7, varests7, numimp, length(unique(dat$Cluster)))
  mmi_upper2[i] <- mmi_uppervec[4]
  
  # Bayesian multilevel multiple imputation
  ests8 <- array(NA, dim = c(numimp, 4))
  varests8 <- array(NA, dim = c(numimp, 4))
  num_clusters <- length(unique(dat2$Cluster))
  size_clusters <- table(dat2$Cluster)
  
  # Initialize and set priors
  numburn <- 1000
  thin <- 100
  numiter <- numburn + thin*numimp
  # Gamma hyperpriors
  c <- d <- 0.01
  
  # Set priors
  # Prior mean for beta
  beta0 <- summary(mmi_impmod)$coef_table[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 12)
  
  # Initial values
  beta <- summary(mmi_impmod)$coef_table[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- dat2$Cluster
  X <- cbind(1, dat2$Y, dat2$A_bin, dat2$X1_bin, dat2$X2_bin,
             dat2$A_bin*dat2$Y,
             dat2$X1_bin*dat2$Y, dat2$A_bin*dat2$X1_bin, dat2$Y*dat2$X2_bin,
             dat2$A_bin*dat2$X2_bin, dat2$A_bin*dat2$Y*dat2$X1_bin,
             dat2$A_bin*dat2$Y*dat2$X2_bin)
  #X <- cbind(1, dat2$Y, dat2$X1_bin, dat2$X2_bin, dat2$X1_bin*dat2$Y,
             #dat2$Y*dat2$X2_bin)
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    datimp8 <- dat2
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(dat2$M)]) / (1 + exp(est[is.na(dat2$M)]))
    datimp8$M[is.na(datimp8$M)] <- rbinom(sum(is.na(dat2$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod_bmmi <- geeglm(Y ~ A_bin*M, family = "gaussian", data = datimp8,
                         id = Cluster, corstr = "exchangeable")
      ests8[(h - numburn) / thin, ] <- coef(mod_bmmi)
      varests8[(h - numburn) / thin, ] <-
        (summary(mod_bmmi)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(dat2)[1], 1, mu)
    z <- (datimp8$M - 1/2) / omega
    
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
  
  bmmi_ests2[i] <- colMeans(ests8)[4]
  bmmi_lowervec <- colMeans(ests8) -
    sqrt(colVars(ests8) + (numimp+1)/(numimp-1)*colMeans(varests8))*
    find_tval(ests8, varests8, numimp, length(unique(dat$Cluster)))
  bmmi_lower2[i] <- bmmi_lowervec[4]
  bmmi_uppervec <- colMeans(ests8) +
    sqrt(colVars(ests8) + (numimp+1)/(numimp-1)*colMeans(varests8))*
    find_tval(ests8, varests8, numimp, length(unique(dat$Cluster)))
  bmmi_upper2[i] <- bmmi_uppervec[4]
  
}

plotdat2 <- data.frame(type = rep("Average"),
                       scenario = rep("Scenario 2", 6),
                       method = c("Real Data", "CCA", "SI", "MI", "MMI",
                                  "B-MMI"),
                       pointest = c(mod$coefficients[4], mean(cca_ests2),
                                    mean(si_ests2), mean(mi_ests2),
                                    mean(mmi_ests2), mean(bmmi_ests2)),
                       lower = c(mod$coefficients[4] -
                                   1.96*summary(mod)$coefficients[4, 2],
                                 mean(cca_lower2),
                                 mean(si_lower2), mean(mi_lower2, na.rm = T),
                                 mean(mmi_lower2, na.rm = T),
                                 mean(bmmi_lower2)),
                       upper = c(mod$coefficients[4] +
                                   1.96*summary(mod)$coefficients[4, 2],
                                 mean(cca_upper2),
                                 mean(si_upper2), mean(mi_upper2, na.rm = T),
                                 mean(mmi_upper2, na.rm = T),
                                 mean(bmmi_upper2)))

plotdat_single2 <- data.frame(type = rep("Single Draw"),
                              scenario = rep("Scenario 2", 6),
                              method = c("Real Data", "CCA", "SI", "MI", "MMI",
                                         "B-MMI"),
                              pointest = c(mod$coefficients[4], cca_ests2[1],
                                           si_ests2[1], mi_ests2[1],
                                           mmi_ests2[1], bmmi_ests2[1]),
                              lower = c(mod$coefficients[4] -
                                          1.96*summary(mod)$coefficients[4, 2],
                                        cca_lower2[1],
                                        si_lower2[1], mi_lower2[1],
                                        mmi_lower2[1],
                                        bmmi_lower2[1]),
                              upper = c(mod$coefficients[4] +
                                          1.96*summary(mod)$coefficients[4, 2],
                                        cca_upper2[1],
                                        si_upper2[1], mi_upper2[1],
                                        mmi_upper2[1],
                                        bmmi_upper2[1]))

ggplot(data = plotdat2, aes(y = method, x = pointest, xmin = lower,
                           xmax = upper)) +
  geom_point() + 
  geom_errorbarh(height = .1) +
  scale_y_discrete(limits = rev(plotdat$method)) +
  labs(x = 'Interaction coefficient', y = 'Method') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5)

ggplot(data = plotdat_single2, aes(y = method, x = pointest, xmin = lower,
                                   xmax = upper)) +
  geom_point() + 
  geom_errorbarh(height = .1) +
  scale_y_discrete(limits = rev(plotdat$method)) +
  labs(x = 'Interaction coefficient', y = 'Method') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5)



#############################################################################
# Scenario 3: Complex MAR with clustered missingness
cca_ests3 <- rep(NA, nsims)
cca_lower3 <- rep(NA, nsims)
cca_upper3 <- rep(NA, nsims)
si_ests3 <- rep(NA, nsims)
si_lower3 <- rep(NA, nsims)
si_upper3 <- rep(NA, nsims)
mi_ests3 <- rep(NA, nsims)
mi_lower3 <- rep(NA, nsims)
mi_upper3 <- rep(NA, nsims)
mmi_ests3 <- rep(NA, nsims)
mmi_lower3 <- rep(NA, nsims)
mmi_upper3 <- rep(NA, nsims)
bmmi_ests3 <- rep(NA, nsims)
bmmi_lower3 <- rep(NA, nsims)
bmmi_upper3 <- rep(NA, nsims)

for (i in 1:nsims) {
  
  dat3 <- dat
  
  ICC_miss <- 0.1
  alpha0_var <- pi^2 * ICC_miss / 3 / (1 - ICC_miss)
  alpha0 <- rnorm(num_clusters, 0, sqrt(alpha0_var))
  
  Rlogit <- 2 + 0.5*dat3$Y - 0.6*dat3$X1_bin - 0.3*dat3$X2_bin +
    0.05*dat3$Y*dat3$X1_bin - 0.15*dat3$Y*dat3$X2_bin +
    #0.1*dat3$Y*dat3$X1_bin*dat3$X2_bin +
    rep(alpha0, size_clusters)
  Rprob <- exp(Rlogit) / (1 + exp(Rlogit))
  dat3$R <- rbinom(dim(dat3)[1], 1, Rprob)
  dat3$M[dat3$R == 0] <- NA
  
  # CCA
  mod_cca <- geeglm(Y ~ A*M, family = "gaussian",
                    data = dat3[complete.cases(dat3), ],
                    id = Cluster, corstr = "exchangeable")
  cca_ests3[i] <- mod_cca$coefficients[4]
  cca_lower3[i] <- mod_cca$coefficients[4] -
    1.96*summary(mod_cca)$coefficients[4, 2]
  cca_upper3[i] <- mod_cca$coefficients[4] +
    1.96*summary(mod_cca)$coefficients[4, 2]
  
  # Single imputation
  impmod <- glm(M ~ Y*A*X1_bin + Y*A*X2_bin, family = "binomial",
                data = dat3[!is.na(dat3$M), ])
  datimp9 <- dat3
  datimp9$M[is.na(datimp9$M)] <- rbinom(length(datimp9$M[is.na(datimp9$M)]), 1,
                                        predict(impmod, dat3[is.na(dat3$M), ],
                                                type = "response"))
  mod_si <- geeglm(Y ~ A*M, family = "gaussian",
                   data = datimp9,
                   id = Cluster, corstr = "exchangeable")
  si_ests3[i] <- mod_si$coefficients[4]
  si_lower3[i] <- mod_si$coefficients[4] -
    1.96*summary(mod_si)$coefficients[4, 2]
  si_upper3[i] <- mod_si$coefficients[4] +
    1.96*summary(mod_si)$coefficients[4, 2]
  
  # Multiple imputation
  numimp <- 15
  ests10 <- array(NA, dim = c(numimp, 4))
  varests10 <- array(NA, dim = c(numimp, 4))
  
  for (m in 1:numimp) {
    
    datimp10 <- dat3
    datimp10$M[is.na(datimp10$M)] <-
      rbinom(length(datimp10$M[is.na(datimp10$M)]), 1,
             predict(impmod, dat3[is.na(dat3$M), ], type = "response"))
    mod_mi <- geeglm(Y ~ A*M, family = "gaussian",
                     data = datimp10,
                     id = Cluster, corstr = "exchangeable")
    ests10[m, ] <- coef(mod_mi)
    varests10[m, ] <- (summary(mod_mi)$coefficients[, 2])^2
    
  }
  
  mi_ests3[i] <- colMeans(ests10)[4]
  mi_lowervec <- colMeans(ests10) -
    sqrt(colVars(ests10) + (numimp+1)/(numimp-1)*colMeans(varests10))*
    find_tval(ests10, varests10, numimp, length(unique(dat$Cluster)))
  mi_lower3[i] <- mi_lowervec[4]
  mi_uppervec <- colMeans(ests10) +
    sqrt(colVars(ests10) + (numimp+1)/(numimp-1)*colMeans(varests10))*
    find_tval(ests10, varests10, numimp, length(unique(dat$Cluster)))
  mi_upper3[i] <- mi_uppervec[4]
  
  # Multilevel multiple imputation
  numimp <- 15
  ests11 <- array(NA, dim = c(numimp, 4))
  varests11 <- array(NA, dim = c(numimp, 4))
  
  mmi_impmod <- mixed_model(M ~ Y*X1_bin + Y*X2_bin, random = ~ 1 | Cluster,
                            family = binomial(), 
                            data = dat2[!is.na(dat2$M), ], iter_EM = 100)
  
  for (m in 1:numimp) {
    
    datimp11 <- dat3
    datimp11$M[is.na(datimp11$M)] <-
      rbinom(length(datimp11$M[is.na(datimp11$M)]), 1,
             predict(mmi_impmod, dat3[is.na(dat3$M), ]))
    mod_mmi <- geeglm(Y ~ A*M, family = "gaussian",
                      data = datimp11,
                      id = Cluster, corstr = "exchangeable")
    ests11[m, ] <- coef(mod_mmi)
    varests11[m, ] <- (summary(mod_mmi)$coefficients[, 2])^2
    
  }
  
  mmi_ests3[i] <- colMeans(ests11)[4]
  mmi_lowervec <- colMeans(ests11) -
    sqrt(colVars(ests11) + (numimp+1)/(numimp-1)*colMeans(varests11))*
    find_tval(ests11, varests11, numimp, length(unique(dat$Cluster)))
  mmi_lower3[i] <- mmi_lowervec[4]
  mmi_uppervec <- colMeans(ests11) +
    sqrt(colVars(ests11) + (numimp+1)/(numimp-1)*colMeans(varests11))*
    find_tval(ests11, varests11, numimp, length(unique(dat$Cluster)))
  mmi_upper3[i] <- mmi_uppervec[4]
  
  # Bayesian multilevel multiple imputation
  ests12 <- array(NA, dim = c(numimp, 4))
  varests12 <- array(NA, dim = c(numimp, 4))
  num_clusters <- length(unique(dat3$Cluster))
  size_clusters <- table(dat3$Cluster)
  
  # Initialize and set priors
  numburn <- 1000
  thin <- 100
  numiter <- numburn + thin*numimp
  # Gamma hyperpriors
  c <- d <- 0.01
  
  # Set priors
  # Prior mean for beta
  beta0 <- summary(mmi_impmod)$coef_table[, 1]
  # Prior precision for beta
  T0 <- diag(.01, 6)
  
  # Initial values
  beta <- summary(mmi_impmod)$coef_table[, 1]
  b <- rep(0, num_clusters) # Random effects
  taub <- 2 # Random effect precision
  
  # Summarize data in helpful vectors and matrices
  #id <- with(df, ave(rep(1, nrow(df)), cluster_ID, FUN = seq_along))
  id <- dat3$Cluster
  X <- cbind(1, dat3$A_bin, dat3$Y, dat3$X1_bin, dat3$X2_bin,
             dat3$A_bin*dat3$Y,
             dat3$X1_bin*dat3$Y, dat3$A_bin*dat3$X1_bin, dat3$Y*dat3$X2_bin,
             dat3$A_bin*dat3$X2_bin, dat3$A_bin*dat3$Y*dat3$X1_bin,
             dat3$A_bin*dat3$Y*dat3$X2_bin)
  X <- cbind(1, dat3$Y, dat3$X1_bin, dat3$X2_bin, dat3$X1_bin*dat3$Y,
             dat3$Y*dat3$X2_bin)
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Impute missing data
    datimp12 <- dat3
    est <- X %*% beta + rep(b, size_clusters)
    estprob <- exp(est[is.na(dat3$M)]) / (1 + exp(est[is.na(dat3$M)]))
    datimp12$M[is.na(datimp12$M)] <- rbinom(sum(is.na(dat3$M)), 1, estprob)
    
    # Fit outcome model after burn-in and after each thinning
    if (h > numburn & h %% thin == 0) {
      mod_bmmi <- geeglm(Y ~ A_bin*M, family = "gaussian", data = datimp12,
                         id = Cluster, corstr = "exchangeable")
      ests12[(h - numburn) / thin, ] <- coef(mod_bmmi)
      varests12[(h - numburn) / thin, ] <-
        (summary(mod_bmmi)$coefficients[, 2])^2
    }
    
    # Update z
    mu <- X %*% beta + rep(b, size_clusters)
    omega <- rpg(dim(dat3)[1], 1, mu)
    z <- (datimp12$M - 1/2) / omega
    
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
  
  bmmi_ests3[i] <- colMeans(ests12)[4]
  bmmi_lowervec <- colMeans(ests12) -
    sqrt(colVars(ests12) + (numimp+1)/(numimp-1)*colMeans(varests12))*
    find_tval(ests12, varests12, numimp, length(unique(dat$Cluster)))
  bmmi_lower3[i] <- bmmi_lowervec[4]
  bmmi_uppervec <- colMeans(ests12) +
    sqrt(colVars(ests12) + (numimp+1)/(numimp-1)*colMeans(varests12))*
    find_tval(ests12, varests12, numimp, length(unique(dat$Cluster)))
  bmmi_upper3[i] <- bmmi_uppervec[4]
  
}

plotdat3 <- data.frame(type = rep("Average"),
                       scenario = rep("Scenario 3", 6),
                       method = c("Real Data", "CCA", "SI", "MI", "MMI",
                                  "B-MMI"),
                       pointest = c(mod$coefficients[4], mean(cca_ests3),
                                    mean(si_ests3), mean(mi_ests3),
                                    mean(mmi_ests3), mean(bmmi_ests3)),
                       lower = c(mod$coefficients[4] -
                                   1.96*summary(mod)$coefficients[4, 2],
                                 mean(cca_lower3),
                                 mean(si_lower3), mean(mi_lower3),
                                 mean(mmi_lower3),
                                 mean(bmmi_lower3)),
                       upper = c(mod$coefficients[4] +
                                   1.96*summary(mod)$coefficients[4, 2],
                                 mean(cca_upper3),
                                 mean(si_upper3), mean(mi_upper3),
                                 mean(mmi_upper3),
                                 mean(bmmi_upper3)))

plotdat_single3 <- data.frame(type = rep("Single Draw"),
                              scenario = rep("Scenario 3", 6),
                              method = c("Real Data", "CCA", "SI", "MI", "MMI",
                                         "B-MMI"),
                              pointest = c(mod$coefficients[4], cca_ests3[1],
                                           si_ests3[1], mi_ests3[1],
                                           mmi_ests3[1], bmmi_ests3[1]),
                              lower = c(mod$coefficients[4] -
                                          1.96*summary(mod)$coefficients[4, 2],
                                        cca_lower3[1],
                                        si_lower3[1], mi_lower3[1],
                                        mmi_lower3[1],
                                        bmmi_lower3[1]),
                              upper = c(mod$coefficients[4] +
                                          1.96*summary(mod)$coefficients[4, 2],
                                        cca_upper3[1],
                                        si_upper3[1], mi_upper3[1],
                                        mmi_upper3[1],
                                        bmmi_upper3[1]))

ggplot(data = plotdat3, aes(y = method, x = pointest, xmin = lower,
                            xmax = upper)) +
  geom_point() + 
  geom_errorbarh(height = .1) +
  scale_y_discrete(limits = rev(plotdat$method)) +
  labs(x = 'Interaction coefficient', y = 'Method') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5)

ggplot(data = plotdat_single3, aes(y = method, x = pointest, xmin = lower,
                                   xmax = upper)) +
  geom_point() + 
  geom_errorbarh(height = .1) +
  scale_y_discrete(limits = rev(plotdat$method)) +
  labs(x = 'Interaction coefficient', y = 'Method') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5)



# Multipanel figure
plotdat_full <- rbind(plotdat_single, plotdat_single2, plotdat_single3,
                      plotdat, plotdat2, plotdat3)

ggplot(data = plotdat_full, aes(y = method, x = pointest, xmin = lower,
                            xmax = upper)) +
  facet_grid(type~scenario) +
  geom_point() + 
  geom_errorbarh(height = .1) +
  scale_y_discrete(limits = rev(plotdat$method)) +
  labs(x = 'Interaction coefficient', y = 'Method') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5)
