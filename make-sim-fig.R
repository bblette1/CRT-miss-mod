rm(list = ls())
library(ggplot2)

# True values of parameters
#true_int <- -0.75
#true_ate <- 1.5 - 0.75*exp(0.5) / (1 + exp(0.5))
pl <- 4 # length of parameter vector

# Make dataset to create figures
makeSimFig <- function(nsims, beta3name, type) {
  
  if (beta3name == 0) {
    true_int <- 0
    true_ate <- 1
  }
  if (beta3name == -1.61) {
    true_int <- -(1 + exp(0.5)) / exp(0.5)
    true_ate <- 0
  }
  
  ###########################################################
  # 20 clusters
  num_clusters <- 20
  results <- array(NA, dim = c(nsims, 228))
  for (i in 1:nsims) {
    if( file.exists(paste("./Results/results", type, "_nc_", num_clusters,
                          "_beta3_", beta3name, "_", i, ".Rdata", sep = ""))) {
      load(paste("./Results/results", type, "_nc_", num_clusters, "_beta3_",
                 beta3name, "_", i, ".Rdata", sep = ""))
      results[i, ] <- sim
    }
  }
  
  # Bias for interaction term
  bias20_int <- colMeans(results[, seq(4, 84, by = 4)], na.rm = T) - true_int
  bias20_int_mcse <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T) /
           (nsims - 1))
  
  # Bias for ATE
  est20_ate <- results[, seq(2, 82, by = 4)]
  bias20_ate <- colMeans(est20_ate, na.rm = T) - true_ate
  bias20_ate_mcse <- sqrt(colMeans((est20_ate - true_ate)^2, na.rm = T) /
                            (nsims - 1))
  
  # Coverage for interaction term
  cov20_int <- rep(NA, 21)
  for (i in 1:6) {
    cov20_int[i] <-
      mean(1*(results[, i*pl] - 1.96*results[, (i + 21)*pl] < true_int &
              results[, i*pl] + 1.96*results[, (i + 21)*pl] > true_int),
           na.rm = T)
  }
  for (i in 7:21) {
    cov20_int[i] <-
      mean(1*(results[, i*pl] -
              results[, (i + 36)*pl]*results[, (i + 21)*pl] < true_int &
              results[, i*pl] +
              results[, (i + 36)*pl]*results[, (i + 21)*pl] > true_int),
           na.rm = T)
  }
  cov20_int_mcse <- sqrt(cov20_int*(1 - cov20_int) / nsims)
  
  # Coverage for ATE
  cov20_ate <- rep(NA, 21)
  for (i in 1:6) {
    cov20_ate[i] <-
      mean(1*(est20_ate[, i] - 1.96*results[, ((i + 21)*pl - 2)] < true_ate &
                est20_ate[, i] + 1.96*results[, ((i + 21)*pl - 2)] > true_ate),
           na.rm = T)
  }
  for (i in 7:21) {
    cov20_ate[i] <-
      mean(1*(est20_ate[, i] - results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)] < true_ate &
                est20_ate[, i] + results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)] > true_ate),
           na.rm = T)
  }
  cov20_ate_mcse <- sqrt(cov20_ate*(1 - cov20_ate) / nsims)
  
  # MSE for interaction term
  mse20_int <-
    colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T)
  mse20_int_mcse <-
    sqrt(colMeans(((results[, seq(4, 84, by = 4)] - true_int)^2 -
                rep(mse20_int, each = 1000))^2, na.rm = T) / (nsims - 1))
  
  # MSE for ATE
  mse20_ate <- colMeans((est20_ate - true_ate)^2, na.rm = T)
  mse20_ate_mcse <-
    sqrt(colMeans(((est20_ate - true_ate)^2 -
                rep(mse20_ate, each = 1000))^2, na.rm = T) / (nsims - 1))
  
  # Power for interaction term
  power20_int <- rep(NA, 21)
  for (i in 1:6) {
    power20_int[i] <-
      mean(1*(results[, i*pl] - 1.96*results[, (i + 21)*pl] > 0 |
                results[, i*pl] + 1.96*results[, (i + 21)*pl] < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    power20_int[i] <-
      mean(1*(results[, i*pl] -
                results[, (i + 36)*pl]*results[, (i + 21)*pl] > 0 |
                results[, i*pl] +
                results[, (i + 36)*pl]*results[, (i + 21)*pl] < 0),
           na.rm = T)
  }
  power20_int_mcse <- sqrt(power20_int*(1 - power20_int) / nsims)
  
  # Power for ATE
  power20_ate <- rep(NA, 21)
  for (i in 1:6) {
    power20_ate[i] <-
      mean(1*(est20_ate[, i] - 1.96*results[, ((i + 21)*pl - 2)] > 0 |
                est20_ate[, i] + 1.96*results[, ((i + 21)*pl - 2)] < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    power20_ate[i] <-
      mean(1*(est20_ate[, i] - results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)]> 0 |
                est20_ate[, i] + results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)]< 0),
           na.rm = T)
  }
  power20_ate_mcse <- sqrt(power20_ate*(1 - power20_ate) / nsims)
  
  ##############################################################
  # Now with 50 clusters
  num_clusters <- 50
  results <- array(NA, dim = c(nsims, 228))
  for (i in 1:nsims) {
    if( file.exists(paste("./Results/results", type, "_nc_", num_clusters,
                          "_beta3_", beta3name, "_", i, ".Rdata", sep = ""))) {
      load(paste("./Results/results", type, "_nc_", num_clusters, "_beta3_",
                 beta3name, "_", i, ".Rdata", sep = ""))
      results[i, ] <- sim
    }
  }
  
  # Bias for interaction term
  bias50_int <- colMeans(results[, seq(4, 84, by = 4)], na.rm = T) - true_int
  bias50_int_mcse <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T) /
           (nsims - 1))
  
  # Bias for ATE
  est50_ate <- results[, seq(2, 82, by = 4)]
  bias50_ate <- colMeans(est50_ate, na.rm = T) - true_ate
  bias50_ate_mcse <- sqrt(colMeans((est50_ate - true_ate)^2, na.rm = T) /
                            (nsims - 1))
  
  # Coverage for interaction term
  cov50_int <- rep(NA, 21)
  for (i in 1:6) {
    cov50_int[i] <-
      mean(1*(results[, i*pl] - 1.96*results[, (i + 21)*pl] < true_int &
                results[, i*pl] + 1.96*results[, (i + 21)*pl] > true_int),
           na.rm = T)
  }
  for (i in 7:21) {
    cov50_int[i] <-
      mean(1*(results[, i*pl] -
                results[, (i + 36)*pl]*results[, (i + 21)*pl] < true_int &
                results[, i*pl] +
                results[, (i + 36)*pl]*results[, (i + 21)*pl] > true_int),
           na.rm = T)
  }
  cov50_int_mcse <- sqrt(cov50_int*(1 - cov50_int) / nsims)
  
  # Coverage for ATE
  cov50_ate <- rep(NA, 21)
  for (i in 1:6) {
    cov50_ate[i] <-
      mean(1*(est50_ate[, i] - 1.96*results[, ((i + 21)*pl - 2)] < true_ate &
                est50_ate[, i] + 1.96*results[, ((i + 21)*pl - 2)] > true_ate),
           na.rm = T)
  }
  for (i in 7:21) {
    cov50_ate[i] <-
      mean(1*(est50_ate[, i] - results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)] < true_ate &
                est50_ate[, i] + results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)] > true_ate),
           na.rm = T)
  }
  cov50_ate_mcse <- sqrt(cov50_ate*(1 - cov50_ate) / nsims)
  
  # RMSE for interaction term
  mse50_int <-
    colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T)
  mse50_int_mcse <-
    sqrt(colMeans(((results[, seq(4, 84, by = 4)] - true_int)^2 -
                     rep(mse50_int, each = 1000))^2, na.rm = T) / (nsims - 1))
  
  # RMSE for ATE
  mse50_ate <-
    colMeans((est50_ate - true_ate)^2, na.rm = T)
  mse50_ate_mcse <-
    sqrt(colMeans(((est50_ate - true_ate)^2 -
                     rep(mse50_ate, each = 1000))^2, na.rm = T) / (nsims - 1))
  
  # Power for interaction term
  power50_int <- rep(NA, 21)
  for (i in 1:6) {
    power50_int[i] <-
      mean(1*(results[, i*pl] - 1.96*results[, (i + 21)*pl] > 0 |
                results[, i*pl] + 1.96*results[, (i + 21)*pl] < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    power50_int[i] <-
      mean(1*(results[, i*pl] -
                results[, (i + 36)*pl]*results[, (i + 21)*pl] > 0 |
                results[, i*pl] +
                results[, (i + 36)*pl]*results[, (i + 21)*pl] < 0),
           na.rm = T)
  }
  power50_int_mcse <- sqrt(power50_int*(1 - power50_int) / nsims)
  
  # Power for ATE
  power50_ate <- rep(NA, 21)
  for (i in 1:6) {
    power50_ate[i] <-
      mean(1*(est50_ate[, i] - 1.96*results[, ((i + 21)*pl - 2)] > 0 |
                est50_ate[, i] + 1.96*results[, ((i + 21)*pl - 2)] < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    power50_ate[i] <-
      mean(1*(est50_ate[, i] - results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)]> 0 |
                est50_ate[, i] + results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)]< 0),
           na.rm = T)
  }
  power50_ate_mcse <- sqrt(power50_ate*(1 - power50_ate) / nsims)
  
  ###############################################################
  # Now with 100 clusters
  num_clusters <- 100
  results <- array(NA, dim = c(nsims, 228))
  for (i in 1:nsims) {
    if( file.exists(paste("./Results/results", type, "_nc_", num_clusters,
                          "_beta3_", beta3name, "_", i, ".Rdata", sep = ""))) {
      load(paste("./Results/results", type, "_nc_", num_clusters, "_beta3_",
                 beta3name, "_", i, ".Rdata", sep = ""))
      results[i, ] <- sim
    }
  }
  
  # Bias for interaction term
  bias100_int <- colMeans(results[, seq(4, 84, by = 4)], na.rm = T) - true_int
  bias100_int_mcse <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T) /
           (nsims - 1))
  
  # Bias for ATE
  est100_ate <- results[, seq(2, 82, by = 4)]
  bias100_ate <- colMeans(est100_ate, na.rm = T) - true_ate
  bias100_ate_mcse <- sqrt(colMeans((est100_ate - true_ate)^2, na.rm = T) /
                            (nsims - 1))
  
  # Coverage for interaction term
  cov100_int <- rep(NA, 21)
  for (i in 1:6) {
    cov100_int[i] <-
      mean(1*(results[, i*pl] - 1.96*results[, (i + 21)*pl] < true_int &
                results[, i*pl] + 1.96*results[, (i + 21)*pl] > true_int),
           na.rm = T)
  }
  for (i in 7:21) {
    cov100_int[i] <-
      mean(1*(results[, i*pl] -
                results[, (i + 36)*pl]*results[, (i + 21)*pl] < true_int &
                results[, i*pl] +
                results[, (i + 36)*pl]*results[, (i + 21)*pl] > true_int),
           na.rm = T)
  }
  cov100_int_mcse <- sqrt(cov100_int*(1 - cov100_int) / nsims)
  
  # Coverage for ATE
  cov100_ate <- rep(NA, 21)
  for (i in 1:6) {
    cov100_ate[i] <-
      mean(1*(est100_ate[, i] - 1.96*results[, ((i + 21)*pl - 2)] < true_ate &
              est100_ate[, i] + 1.96*results[, ((i + 21)*pl - 2)] > true_ate),
           na.rm = T)
  }
  for (i in 7:21) {
    cov100_ate[i] <-
      mean(1*(est100_ate[, i] - results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)] < true_ate &
                est100_ate[, i] + results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)] > true_ate),
           na.rm = T)
  }
  cov100_ate_mcse <- sqrt(cov100_ate*(1 - cov100_ate) / nsims)
  
  # MSE for interaction term
  mse100_int <-
    colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T)
  mse100_int_mcse <-
    sqrt(colMeans(((results[, seq(4, 84, by = 4)] - true_int)^2 -
                     rep(mse100_int, each = 1000))^2, na.rm = T) / (nsims - 1))
  
  # MSE for ATE
  mse100_ate <-
    colMeans((est100_ate - true_ate)^2, na.rm = T)
  mse100_ate_mcse <-
    sqrt(colMeans(((est100_ate - true_ate)^2 -
                     rep(mse100_ate, each = 1000))^2, na.rm = T) / (nsims - 1))
  
  # Power for interaction term
  power100_int <- rep(NA, 21)
  for (i in 1:6) {
    power100_int[i] <-
      mean(1*(results[, i*pl] - 1.96*results[, (i + 21)*pl] > 0 |
                results[, i*pl] + 1.96*results[, (i + 21)*pl] < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    power100_int[i] <-
      mean(1*(results[, i*pl] -
                results[, (i + 36)*pl]*results[, (i + 21)*pl] > 0 |
                results[, i*pl] +
                results[, (i + 36)*pl]*results[, (i + 21)*pl] < 0),
           na.rm = T)
  }
  power100_int_mcse <- sqrt(power100_int*(1 - power100_int) / nsims)
  
  # Power for ATE
  power100_ate <- rep(NA, 21)
  for (i in 1:6) {
    power100_ate[i] <-
      mean(1*(est100_ate[, i] - 1.96*results[, ((i + 21)*pl - 2)] > 0 |
                est100_ate[, i] + 1.96*results[, ((i + 21)*pl - 2)] < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    power100_ate[i] <-
      mean(1*(est100_ate[, i] - results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)]> 0 |
                est100_ate[, i] + results[, ((i + 36)*pl - 2)]*
                results[, ((i + 21)*pl - 2)]< 0),
           na.rm = T)
  }
  power100_ate_mcse <- sqrt(power100_ate*(1 - power100_ate) / nsims)
  
  # Summary figure
  fig_dat <-
    data.frame("Output" = c(rep(bias20_int[1], 4), bias20_int[-1],
                            rep(bias50_int[1], 4), bias50_int[-1],
                            rep(bias100_int[1], 4), bias100_int[-1],
                            rep(cov20_int[1], 4), cov20_int[-1],
                            rep(cov50_int[1], 4), cov50_int[-1],
                            rep(cov100_int[1], 4), cov100_int[-1],
                            rep(mse20_int[1], 4), mse20_int[-1],
                            rep(mse50_int[1], 4), mse50_int[-1],
                            rep(mse100_int[1], 4), mse100_int[-1],
                            rep(power20_int[1], 4), power20_int[-1],
                            rep(power50_int[1], 4), power50_int[-1],
                            rep(power100_int[1], 4), power100_int[-1],
                            rep(bias20_ate[1], 4), bias20_ate[-1],
                            rep(bias50_ate[1], 4), bias50_ate[-1],
                            rep(bias100_ate[1], 4), bias100_ate[-1],
                            rep(cov20_ate[1], 4), cov20_ate[-1],
                            rep(cov50_ate[1], 4), cov50_ate[-1],
                            rep(cov100_ate[1], 4), cov100_ate[-1],
                            rep(mse20_ate[1], 4), mse20_ate[-1],
                            rep(mse50_ate[1], 4), mse50_ate[-1],
                            rep(mse100_ate[1], 4), mse100_ate[-1],
                            rep(power20_ate[1], 4), power20_ate[-1],
                            rep(power50_ate[1], 4), power50_ate[-1],
                            rep(power100_ate[1], 4), power100_ate[-1]),
               "Estimand" = rep(c("Int", "ATE"), each = 288),
               "Metric" = rep(rep(c("Bias", "Coverage", "MSE", "Power"),
                                  each = 72),
                              2),
               "Hline" = c(rep(c(0, 0.95, 0, (1 - 0.95*(beta3name == 0))),
                               each = 72),
                           rep(c(0, 0.95, 0, (1 - 0.95*(beta3name != 0))),
                               each = 72)),
               "Method" = rep(factor(rep(c("SI", "MI", "MMI", "B-MMI",
                                          rep(c("SI", "MI", "MMI", "B-MMI"),
                                       each = 5)), 12),
                                    levels = c("SI", "MI", 
                                               "MMI", "B-MMI")), 2),
               "Model" = rep(rep(c("CCA", "CCA", "CCA", "CCA",
                               rep(c("M~X+A+Y", "M~X+A*Y", "M~X*A+Y",
                                     "M~X*A+Y*A", "M~X*A*Y"),
                                   4)), 12), 2),
               "Num_Clusters" = rep(rep(rep(c(20, 50, 100), each = 24), 4), 2))
  
  mcse_table <- data.frame("Output" = c(bias20_int_mcse, bias50_int_mcse,
                                        bias100_int_mcse, cov20_int_mcse,
                                        cov50_int_mcse, cov100_int_mcse,
                                        mse20_int_mcse, mse50_int_mcse,
                                        mse100_int_mcse, power20_int_mcse,
                                        power50_int_mcse, power100_int_mcse,
                                        bias20_ate_mcse, bias50_ate_mcse,
                                        bias100_ate_mcse, cov20_ate_mcse,
                                        cov50_ate_mcse, cov100_ate_mcse,
                                        mse20_ate_mcse, mse50_ate_mcse,
                                        mse100_ate_mcse, power20_ate_mcse,
                                        power50_ate_mcse, power100_ate_mcse),
                           "Estimand" = rep(c("Int", "Ate"), each = 12*21),
                           "Metric" = rep(rep(c("Bias", "Coverage", "MSE",
                                            "Power"), each = 3*21), 2),
                           "Method_Model" = c(1))
  
  return(list(fig_dat, mcse_table))
  
}


# Plot first simulation scenario results
figdat1_0 <- makeSimFig(2000, 0, "")[[1]]
figdat1_1.61 <- makeSimFig(2000, -1.61, "")[[1]]

# MCSE look appropriate
#tabdat1_0 <- makeSimFig(2000, 0)[[2]]

# Beta3 = 0
newlabs <- c("Bias", "Coverage", "MSE", "Type I Error")
names(newlabs) <- c("Bias", "Coverage", "MSE", "Power")
ggplot(data = figdat1_0[figdat1_0$Estimand == "Int", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free",
             labeller = labeller(Metric = newlabs)) +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))

ggplot(data = figdat1_0[figdat1_0$Estimand == "ATE", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))

# Beta3 = -1.61
ggplot(data = figdat1_1.61[figdat1_1.61$Estimand == "Int", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))

ggplot(data = figdat1_1.61[figdat1_1.61$Estimand == "ATE", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free",
             labeller = labeller(Metric = newlabs)) +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))



# Plot second simulation scenario results
figdat2_0 <- makeSimFig(2000, 0, "2")[[1]]
figdat2_1.61 <- makeSimFig(2000, -1.61, "2")[[1]]

ggplot(data = figdat2_0[figdat2_0$Estimand == "Int", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free",
             labeller = labeller(Metric = newlabs)) +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))

ggplot(data = figdat2_0[figdat2_0$Estimand == "ATE", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))

# Beta3 = -1.61
ggplot(data = figdat2_1.61[figdat2_1.61$Estimand == "Int", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))

ggplot(data = figdat2_1.61[figdat2_1.61$Estimand == "ATE", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(Metric ~ Method, scales = "free",
             labeller = labeller(Metric = newlabs)) +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(strip.background = element_rect(fill = "dark blue"))
