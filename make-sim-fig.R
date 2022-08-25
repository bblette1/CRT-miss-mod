rm(list = ls())
library(ggplot2)

# True values of parameters
true_int <- -0.75
true_ate <- 1.5 - 0.75*exp(0.5) / (1 + exp(0.5))
pl <- 4 # length of parameter vector

# Make dataset to create figures
makeSimFig <- function(type, nsims, ICC_out, ICC_mod) {
  
  ###########################################################
  # 20 clusters
  num_clusters <- 20
  results <- array(NA, dim = c(nsims, 228))
  for (i in 1:nsims) {
    if( file.exists(paste("./Results/results_", type, "_Iout_", ICC_out,
                          "_Imod_", ICC_mod, "_nc_",
                          num_clusters, "_", i, ".Rdata", sep = ""))) {
      load(paste("./Results/results_", type, "_Iout_", ICC_out,
                 "_Imod_", ICC_mod, "_nc_",
                 num_clusters, "_", i, ".Rdata", sep = ""))
      results[i, ] <- sim
    }
  }
  
  # Bias for interaction term
  bias20_int <- colMeans(results[, seq(4, 84, by = 4)], na.rm = T) - true_int
  
  # Bias for ATE
  est20_ate <- results[, seq(2, 82, by = 4)] +
               exp(0.5) / (1 + exp(0.5)) * results[, seq(4, 84, by = 4)]
  bias20_ate <- colMeans(est20_ate, na.rm = T) - true_ate
  
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
  
  # Coverage for ATE
  cov20_ate <- rep(NA, 21)
  for (i in 1:6) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    cov20_ate[i] <-
      mean(1*(est20_ate[, i] - 1.96*se < true_ate &
                est20_ate[, i] + 1.96*se > true_ate),
           na.rm = T)
  }
  for (i in 7:21) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    cov20_ate[i] <-
      mean(1*(est20_ate[, i] - results[, (i + 36)*pl]*se < true_ate &
                est20_ate[, i] + results[, (i + 36)*pl]*se > true_ate),
           na.rm = T)
  }
  
  # RMSE for interaction term
  rmse20_int <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T))
  
  # RMSE for ATE
  rmse20_ate <-
    sqrt(colMeans((est20_ate - true_ate)^2, na.rm = T))
  
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
  
  # Power for ATE
  power20_ate <- rep(NA, 21)
  for (i in 1:6) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    power20_ate[i] <-
      mean(1*(est20_ate[, i] - 1.96*se > 0 |
                est20_ate[, i] + 1.96*se < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    power20_ate[i] <-
      mean(1*(est20_ate[, i] - results[, (i + 36)*pl]*se > 0 |
                est20_ate[, i] + results[, (i + 36)*pl]*se < 0),
           na.rm = T)
  }
  
  ##############################################################
  # Now with 50 clusters
  num_clusters <- 50
  results <- array(NA, dim = c(nsims, 228))
  for (i in 1:nsims) {
    if( file.exists(paste("./Results/results_", type, "_Iout_", ICC_out,
                          "_Imod_", ICC_mod, "_nc_",
                          num_clusters, "_", i, ".Rdata", sep = ""))) {
      load(paste("./Results/results_", type, "_Iout_", ICC_out,
                 "_Imod_", ICC_mod, "_nc_",
                 num_clusters, "_", i, ".Rdata", sep = ""))
      results[i, ] <- sim
    }
  }
  
  # Bias for interaction term
  bias50_int <- colMeans(results[, seq(4, 84, by = 4)], na.rm = T) - true_int
  
  # Bias for ATE
  est50_ate <- results[, seq(2, 82, by = 4)] +
               exp(0.5) / (1 + exp(0.5)) * results[, seq(4, 84, by = 4)]
  bias50_ate <- colMeans(est50_ate, na.rm = T) - true_ate
  
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
  
  # Coverage for ATE
  cov50_ate <- rep(NA, 21)
  for (i in 1:6) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    cov50_ate[i] <-
      mean(1*(est50_ate[, i] - 1.96*se < true_ate &
                est50_ate[, i] + 1.96*se > true_ate),
           na.rm = T)
  }
  for (i in 7:21) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    cov50_ate[i] <-
      mean(1*(est50_ate[, i] - results[, (i + 36)*pl]*se < true_ate &
                est50_ate[, i] + results[, (i + 36)*pl]*se > true_ate),
           na.rm = T)
  }
  
  # RMSE for interaction term
  rmse50_int <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T))
  
  # RMSE for ATE
  rmse50_ate <-
    sqrt(colMeans((est50_ate - true_ate)^2, na.rm = T))
  
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
  
  # Power for ATE
  power50_ate <- rep(NA, 21)
  for (i in 1:6) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    power50_ate[i] <-
      mean(1*(est50_ate[, i] - 1.96*se > 0 |
                est50_ate[, i] + 1.96*se < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    power50_ate[i] <-
      mean(1*(est50_ate[, i] - results[, (i + 36)*pl]*se > 0 |
                est50_ate[, i] + results[, (i + 36)*pl]*se < 0),
           na.rm = T)
  }
  
  ###############################################################
  # Now with 100 clusters
  num_clusters <- 100
  results <- array(NA, dim = c(nsims, 228))
  for (i in 1:nsims) {
    if( file.exists(paste("./Results/results_", type, "_Iout_", ICC_out,
                          "_Imod_", ICC_mod, "_nc_",
                          num_clusters, "_", i, ".Rdata", sep = ""))) {
      load(paste("./Results/results_", type, "_Iout_", ICC_out,
                 "_Imod_", ICC_mod, "_nc_",
                 num_clusters, "_", i, ".Rdata", sep = ""))
      results[i, ] <- sim
    }
  }
  
  # Bias for interaction term
  bias100_int <- colMeans(results[, seq(4, 84, by = 4)], na.rm = T) - true_int
  
  # Bias for ATE
  est100_ate <- results[, seq(2, 82, by = 4)] +
                exp(0.5) / (1 + exp(0.5)) * results[, seq(4, 84, by = 4)]
  bias100_ate <- colMeans(est100_ate, na.rm = T) - true_ate
  
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
  
  # Coverage for ATE
  cov100_ate <- rep(NA, 21)
  for (i in 1:6) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    cov100_ate[i] <-
      mean(1*(est100_ate[, i] - 1.96*se < true_ate &
                est100_ate[, i] + 1.96*se > true_ate),
           na.rm = T)
  }
  for (i in 7:21) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    cov100_ate[i] <-
      mean(1*(est100_ate[, i] - results[, (i + 36)*pl]*se < true_ate &
                est100_ate[, i] + results[, (i + 36)*pl]*se > true_ate),
           na.rm = T)
  }
  
  # RMSE for interaction term
  rmse100_int <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T))
  
  # RMSE for ATE
  rmse100_ate <-
    sqrt(colMeans((est100_ate - true_ate)^2, na.rm = T))
  
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
  
  # Power for ATE
  power100_ate <- rep(NA, 21)
  for (i in 1:6) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    power100_ate[i] <-
      mean(1*(est100_ate[, i] - 1.96*se > 0 |
                est100_ate[, i] + 1.96*se < 0),
           na.rm = T)
  }
  for (i in 7:21) {
    se <- sqrt((results[, ((i + 21)*pl - 2)])^2 +
                 (0.622^2)*(results[, (i + 21)*pl])^2)
    power100_ate[i] <-
      mean(1*(est100_ate[, i] - results[, (i + 36)*pl]*se > 0 |
                est100_ate[, i] + results[, (i + 36)*pl]*se < 0),
           na.rm = T)
  }
  
  # Summary figure
  fig_dat_int <-
    data.frame("Output" = c(rep(bias20_int[1], 4), bias20_int[-1],
                            rep(bias50_int[1], 4), bias50_int[-1],
                            rep(bias100_int[1], 4), bias100_int[-1],
                            rep(cov20_int[1], 4), cov20_int[-1],
                            rep(cov50_int[1], 4), cov50_int[-1],
                            rep(cov100_int[1], 4), cov100_int[-1],
                            rep(rmse20_int[1], 4), rmse20_int[-1],
                            rep(rmse50_int[1], 4), rmse50_int[-1],
                            rep(rmse100_int[1], 4), rmse100_int[-1],
                            rep(power20_int[1], 4), power20_int[-1],
                            rep(power50_int[1], 4), power50_int[-1],
                            rep(power100_int[1], 4), power100_int[-1],
                            rep(bias20_ate[1], 4), bias20_ate[-1],
                            rep(bias50_ate[1], 4), bias50_ate[-1],
                            rep(bias100_ate[1], 4), bias100_ate[-1],
                            rep(cov20_ate[1], 4), cov20_ate[-1],
                            rep(cov50_ate[1], 4), cov50_ate[-1],
                            rep(cov100_ate[1], 4), cov100_ate[-1],
                            rep(rmse20_ate[1], 4), rmse20_ate[-1],
                            rep(rmse50_ate[1], 4), rmse50_ate[-1],
                            rep(rmse100_ate[1], 4), rmse100_ate[-1],
                            rep(power20_ate[1], 4), power20_ate[-1],
                            rep(power50_ate[1], 4), power50_ate[-1],
                            rep(power100_ate[1], 4), power100_ate[-1]),
               "Estimand" = rep(c("Int", "ATE"), each = 288),
               "Metric" = rep(rep(c("Bias", "Coverage", "RMSE", "Power"),
                                  each = 72),
                              2),
               "Hline" = rep(rep(c(0, 0.95, 0, 0), each = 72), 2),
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
  
  return(fig_dat_int)
  
}


# Plot first simulation scenario results
figdat1 <- makeSimFig("out", 1000, 0.1, 0.1)

ggplot(data = figdat1[figdat1$Estimand == "Int", ],
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

ggplot(data = figdat1[figdat1$Estimand == "ATE", ],
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



# Plot second simulation scenario results
figdat2 <- makeSimFig("out2", 1000, 0.1, 0.1)

ggplot(data = figdat2[figdat2$Estimand == "Int", ],
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

ggplot(data = figdat2[figdat2$Estimand == "ATE", ],
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
