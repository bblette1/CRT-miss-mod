# True values of parameters
true_int <- -0.75
true_ate <- 1.5 - 0.75*exp(0.5) / (1 + exp(0.5))
pl <- 4 # length of parameter vector

# Make dataset to create figures
makeSimFig <- function(type, nsims, ICC_out, ICC_mod) {
  
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
  
  # Bias for adjusted treatment effect
  est20_ate <- colMeans(results[, seq(2, 82, by = 4)], na.rm = T) +
    exp(0.5) / (1 + exp(0.5))*
    colMeans(results[, seq(4, 84, by = 4)], na.rm = T)
  bias20_ate <- est20_ate - true_ate
  
  # Coverage for interaction term
  cov20_int_1 <- 1*(results[, 1*pl] - 1.96*results[, 22*pl] < true_int &
                      results[, 1*pl] + 1.96*results[, 22*pl] > true_int)
  cov20_int_2 <- 1*(results[, 2*pl] - 1.96*results[, 23*pl] < true_int &
                      results[, 2*pl] + 1.96*results[, 23*pl] > true_int)
  cov20_int_3 <- 1*(results[, 3*pl] - 1.96*results[, 24*pl] < true_int &
                      results[, 3*pl] + 1.96*results[, 24*pl] > true_int)
  cov20_int_4 <- 1*(results[, 4*pl] - 1.96*results[, 25*pl] < true_int &
                      results[, 4*pl] + 1.96*results[, 25*pl] > true_int)
  cov20_int_5 <- 1*(results[, 5*pl] - 1.96*results[, 26*pl] < true_int &
                      results[, 5*pl] + 1.96*results[, 26*pl] > true_int)
  cov20_int_6 <- 1*(results[, 6*pl] - 1.96*results[, 27*pl] < true_int &
                      results[, 6*pl] + 1.96*results[, 27*pl] > true_int)
  cov20_int_7 <-
    1*(results[, 7*pl] - results[, 43*pl]*results[, 28*pl] < true_int &
         results[, 7*pl] + results[, 43*pl]*results[, 28*pl] > true_int)
  cov20_int_8 <-
    1*(results[, 8*pl] - results[, 44*pl]*results[, 29*pl] < true_int &
         results[, 8*pl] + results[, 44*pl]*results[, 29*pl] > true_int)
  cov20_int_9 <-
    1*(results[, 9*pl] - results[, 45*pl]*results[, 30*pl] < true_int &
         results[, 9*pl] + results[, 45*pl]*results[, 30*pl] > true_int)
  cov20_int_10 <-
    1*(results[, 10*pl] - results[, 46*pl]*results[, 31*pl] < true_int &
         results[, 10*pl] + results[, 46*pl]*results[, 31*pl] > true_int)
  cov20_int_11 <-
    1*(results[, 11*pl] - results[, 47*pl]*results[, 32*pl] < true_int &
         results[, 11*pl] + results[, 47*pl]*results[, 32*pl] > true_int)
  cov20_int_12 <-
    1*(results[, 12*pl] - results[, 48*pl]*results[, 33*pl] < true_int &
         results[, 12*pl] + results[, 48*pl]*results[, 33*pl] > true_int)
  cov20_int_13 <-
    1*(results[, 13*pl] - results[, 49*pl]*results[, 34*pl] < true_int &
         results[, 13*pl] + results[, 49*pl]*results[, 34*pl] > true_int)
  cov20_int_14 <-
    1*(results[, 14*pl] - results[, 50*pl]*results[, 35*pl] < true_int &
         results[, 14*pl] + results[, 50*pl]*results[, 35*pl] > true_int)
  cov20_int_15 <-
    1*(results[, 15*pl] - results[, 51*pl]*results[, 36*pl] < true_int &
         results[, 15*pl] + results[, 51*pl]*results[, 36*pl] > true_int)
  cov20_int_16 <-
    1*(results[, 16*pl] - results[, 52*pl]*results[, 37*pl] < true_int &
         results[, 16*pl] + results[, 52*pl]*results[, 37*pl] > true_int)
  cov20_int_17 <-
    1*(results[, 17*pl] - results[, 53*pl]*results[, 38*pl] < true_int &
         results[, 17*pl] + results[, 53*pl]*results[, 38*pl] > true_int)
  cov20_int_18 <-
    1*(results[, 18*pl] - results[, 54*pl]*results[, 39*pl] < true_int &
         results[, 18*pl] + results[, 54*pl]*results[, 39*pl] > true_int)
  cov20_int_19 <-
    1*(results[, 19*pl] - results[, 55*pl]*results[, 40*pl] < true_int &
         results[, 19*pl] + results[, 55*pl]*results[, 40*pl] > true_int)
  cov20_int_20 <-
    1*(results[, 20*pl] - results[, 56*pl]*results[, 41*pl] < true_int &
         results[, 20*pl] + results[, 56*pl]*results[, 41*pl] > true_int)
  cov20_int_21 <-
    1*(results[, 21*pl] - results[, 57*pl]*results[, 42*pl] < true_int &
         results[, 21*pl] + results[, 57*pl]*results[, 42*pl] > true_int)
  
  # Coverage for adjusted treatment effect
  cov20_ate_1 <- 1*(est20_ate[1] - 1.96*results[, 22*pl] < true_int &
                      results[, 1*pl] + 1.96*results[, 22*pl] > true_int)
  
  # RMSE for interaction term
  rmse20_int <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T))
  
  # RMSE for adjusted treatment effect
  
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
  
  # Bias for adjusted treatment effect
  est50_ate <- colMeans(results[, seq(2, 82, by = 4)], na.rm = T) +
    exp(0.5) / (1 + exp(0.5))*
    colMeans(results[, seq(4, 84, by = 4)], na.rm = T)
  bias50_ate <- est50_ate - true_ate
  
  # Coverage for interaction term
  cov50_int_1 <- 1*(results[, 1*pl] - 1.96*results[, 22*pl] < true_int &
                      results[, 1*pl] + 1.96*results[, 22*pl] > true_int)
  cov50_int_2 <- 1*(results[, 2*pl] - 1.96*results[, 23*pl] < true_int &
                      results[, 2*pl] + 1.96*results[, 23*pl] > true_int)
  cov50_int_3 <- 1*(results[, 3*pl] - 1.96*results[, 24*pl] < true_int &
                      results[, 3*pl] + 1.96*results[, 24*pl] > true_int)
  cov50_int_4 <- 1*(results[, 4*pl] - 1.96*results[, 25*pl] < true_int &
                      results[, 4*pl] + 1.96*results[, 25*pl] > true_int)
  cov50_int_5 <- 1*(results[, 5*pl] - 1.96*results[, 26*pl] < true_int &
                      results[, 5*pl] + 1.96*results[, 26*pl] > true_int)
  cov50_int_6 <- 1*(results[, 6*pl] - 1.96*results[, 27*pl] < true_int &
                      results[, 6*pl] + 1.96*results[, 27*pl] > true_int)
  cov50_int_7 <-
    1*(results[, 7*pl] - results[, 43*pl]*results[, 28*pl] < true_int &
         results[, 7*pl] + results[, 43*pl]*results[, 28*pl] > true_int)
  cov50_int_8 <-
    1*(results[, 8*pl] - results[, 44*pl]*results[, 29*pl] < true_int &
         results[, 8*pl] + results[, 44*pl]*results[, 29*pl] > true_int)
  cov50_int_9 <-
    1*(results[, 9*pl] - results[, 45*pl]*results[, 30*pl] < true_int &
         results[, 9*pl] + results[, 45*pl]*results[, 30*pl] > true_int)
  cov50_int_10 <-
    1*(results[, 10*pl] - results[, 46*pl]*results[, 31*pl] < true_int &
         results[, 10*pl] + results[, 46*pl]*results[, 31*pl] > true_int)
  cov50_int_11 <-
    1*(results[, 11*pl] - results[, 47*pl]*results[, 32*pl] < true_int &
         results[, 11*pl] + results[, 47*pl]*results[, 32*pl] > true_int)
  cov50_int_12 <-
    1*(results[, 12*pl] - results[, 48*pl]*results[, 33*pl] < true_int &
         results[, 12*pl] + results[, 48*pl]*results[, 33*pl] > true_int)
  cov50_int_13 <-
    1*(results[, 13*pl] - results[, 49*pl]*results[, 34*pl] < true_int &
         results[, 13*pl] + results[, 49*pl]*results[, 34*pl] > true_int)
  cov50_int_14 <-
    1*(results[, 14*pl] - results[, 50*pl]*results[, 35*pl] < true_int &
         results[, 14*pl] + results[, 50*pl]*results[, 35*pl] > true_int)
  cov50_int_15 <-
    1*(results[, 15*pl] - results[, 51*pl]*results[, 36*pl] < true_int &
         results[, 15*pl] + results[, 51*pl]*results[, 36*pl] > true_int)
  cov50_int_16 <-
    1*(results[, 16*pl] - results[, 52*pl]*results[, 37*pl] < true_int &
         results[, 16*pl] + results[, 52*pl]*results[, 37*pl] > true_int)
  cov50_int_17 <-
    1*(results[, 17*pl] - results[, 53*pl]*results[, 38*pl] < true_int &
         results[, 17*pl] + results[, 53*pl]*results[, 38*pl] > true_int)
  cov50_int_18 <-
    1*(results[, 18*pl] - results[, 54*pl]*results[, 39*pl] < true_int &
         results[, 18*pl] + results[, 54*pl]*results[, 39*pl] > true_int)
  cov50_int_19 <-
    1*(results[, 19*pl] - results[, 55*pl]*results[, 40*pl] < true_int &
         results[, 19*pl] + results[, 55*pl]*results[, 40*pl] > true_int)
  cov50_int_20 <-
    1*(results[, 20*pl] - results[, 56*pl]*results[, 41*pl] < true_int &
         results[, 20*pl] + results[, 56*pl]*results[, 41*pl] > true_int)
  cov50_int_21 <-
    1*(results[, 21*pl] - results[, 57*pl]*results[, 42*pl] < true_int &
         results[, 21*pl] + results[, 57*pl]*results[, 42*pl] > true_int)
  
  # Coverage for adjusted treatment effect
  cov50_ate_1 <- 1*(est50_ate[1] - 1.96*results[, 22*pl] < true_int &
                      results[, 1*pl] + 1.96*results[, 22*pl] > true_int)
  
  # RMSE for interaction term
  rmse50_int <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T))
  
  # RMSE for adjusted treatment effect

  
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
  
  # Bias for adjusted treatment effect
  est100_ate <- colMeans(results[, seq(2, 82, by = 4)], na.rm = T) +
    exp(0.5) / (1 + exp(0.5))*
    colMeans(results[, seq(4, 84, by = 4)], na.rm = T)
  bias100_ate <- est100_ate - true_ate
  
  # Coverage for interaction term
  cov100_int_1 <- 1*(results[, 1*pl] - 1.96*results[, 22*pl] < true_int &
                      results[, 1*pl] + 1.96*results[, 22*pl] > true_int)
  cov100_int_2 <- 1*(results[, 2*pl] - 1.96*results[, 23*pl] < true_int &
                      results[, 2*pl] + 1.96*results[, 23*pl] > true_int)
  cov100_int_3 <- 1*(results[, 3*pl] - 1.96*results[, 24*pl] < true_int &
                      results[, 3*pl] + 1.96*results[, 24*pl] > true_int)
  cov100_int_4 <- 1*(results[, 4*pl] - 1.96*results[, 25*pl] < true_int &
                      results[, 4*pl] + 1.96*results[, 25*pl] > true_int)
  cov100_int_5 <- 1*(results[, 5*pl] - 1.96*results[, 26*pl] < true_int &
                      results[, 5*pl] + 1.96*results[, 26*pl] > true_int)
  cov100_int_6 <- 1*(results[, 6*pl] - 1.96*results[, 27*pl] < true_int &
                      results[, 6*pl] + 1.96*results[, 27*pl] > true_int)
  cov100_int_7 <-
    1*(results[, 7*pl] - results[, 43*pl]*results[, 28*pl] < true_int &
         results[, 7*pl] + results[, 43*pl]*results[, 28*pl] > true_int)
  cov100_int_8 <-
    1*(results[, 8*pl] - results[, 44*pl]*results[, 29*pl] < true_int &
         results[, 8*pl] + results[, 44*pl]*results[, 29*pl] > true_int)
  cov100_int_9 <-
    1*(results[, 9*pl] - results[, 45*pl]*results[, 30*pl] < true_int &
         results[, 9*pl] + results[, 45*pl]*results[, 30*pl] > true_int)
  cov100_int_10 <-
    1*(results[, 10*pl] - results[, 46*pl]*results[, 31*pl] < true_int &
         results[, 10*pl] + results[, 46*pl]*results[, 31*pl] > true_int)
  cov100_int_11 <-
    1*(results[, 11*pl] - results[, 47*pl]*results[, 32*pl] < true_int &
         results[, 11*pl] + results[, 47*pl]*results[, 32*pl] > true_int)
  cov100_int_12 <-
    1*(results[, 12*pl] - results[, 48*pl]*results[, 33*pl] < true_int &
         results[, 12*pl] + results[, 48*pl]*results[, 33*pl] > true_int)
  cov100_int_13 <-
    1*(results[, 13*pl] - results[, 49*pl]*results[, 34*pl] < true_int &
         results[, 13*pl] + results[, 49*pl]*results[, 34*pl] > true_int)
  cov100_int_14 <-
    1*(results[, 14*pl] - results[, 50*pl]*results[, 35*pl] < true_int &
         results[, 14*pl] + results[, 50*pl]*results[, 35*pl] > true_int)
  cov100_int_15 <-
    1*(results[, 15*pl] - results[, 51*pl]*results[, 36*pl] < true_int &
         results[, 15*pl] + results[, 51*pl]*results[, 36*pl] > true_int)
  cov100_int_16 <-
    1*(results[, 16*pl] - results[, 52*pl]*results[, 37*pl] < true_int &
         results[, 16*pl] + results[, 52*pl]*results[, 37*pl] > true_int)
  cov100_int_17 <-
    1*(results[, 17*pl] - results[, 53*pl]*results[, 38*pl] < true_int &
         results[, 17*pl] + results[, 53*pl]*results[, 38*pl] > true_int)
  cov100_int_18 <-
    1*(results[, 18*pl] - results[, 54*pl]*results[, 39*pl] < true_int &
         results[, 18*pl] + results[, 54*pl]*results[, 39*pl] > true_int)
  cov100_int_19 <-
    1*(results[, 19*pl] - results[, 55*pl]*results[, 40*pl] < true_int &
         results[, 19*pl] + results[, 55*pl]*results[, 40*pl] > true_int)
  cov100_int_20 <-
    1*(results[, 20*pl] - results[, 56*pl]*results[, 41*pl] < true_int &
         results[, 20*pl] + results[, 56*pl]*results[, 41*pl] > true_int)
  cov100_int_21 <-
    1*(results[, 21*pl] - results[, 57*pl]*results[, 42*pl] < true_int &
         results[, 21*pl] + results[, 57*pl]*results[, 42*pl] > true_int)
  
  # Coverage for adjusted treatment effect
  cov100_ate_1 <- 1*(est100_ate[1] - 1.96*results[, 22*pl] < true_int &
                      results[, 1*pl] + 1.96*results[, 22*pl] > true_int)
  
  # RMSE for interaction term
  rmse100_int <-
    sqrt(colMeans((results[, seq(4, 84, by = 4)] - true_int)^2, na.rm = T))
  
  # RMSE for adjusted treatment effect
  
  
  # Summary figure
  fig_dat_int <-
    data.frame("Output" = c(bias20_int, bias50_int, bias100_int,
                            mean(cov20_int_1, na.rm = T),
                            mean(cov20_int_2, na.rm = T),
                            mean(cov20_int_3, na.rm = T),
                            mean(cov20_int_4, na.rm = T),
                            mean(cov20_int_5, na.rm = T),
                            mean(cov20_int_6, na.rm = T),
                            mean(cov20_int_7, na.rm = T),
                            mean(cov20_int_8, na.rm = T),
                            mean(cov20_int_9, na.rm = T),
                            mean(cov20_int_10, na.rm = T),
                            mean(cov20_int_11, na.rm = T),
                            mean(cov20_int_12, na.rm = T),
                            mean(cov20_int_13, na.rm = T),
                            mean(cov20_int_14, na.rm = T),
                            mean(cov20_int_15, na.rm = T),
                            mean(cov20_int_16, na.rm = T),
                            mean(cov20_int_17, na.rm = T),
                            mean(cov20_int_18, na.rm = T),
                            mean(cov20_int_19, na.rm = T),
                            mean(cov20_int_20, na.rm = T),
                            mean(cov20_int_21, na.rm = T),
                            mean(cov50_int_1, na.rm = T),
                            mean(cov50_int_2, na.rm = T),
                            mean(cov50_int_3, na.rm = T),
                            mean(cov50_int_4, na.rm = T),
                            mean(cov50_int_5, na.rm = T),
                            mean(cov50_int_6, na.rm = T),
                            mean(cov50_int_7, na.rm = T),
                            mean(cov50_int_8, na.rm = T),
                            mean(cov50_int_9, na.rm = T),
                            mean(cov50_int_10, na.rm = T),
                            mean(cov50_int_11, na.rm = T),
                            mean(cov50_int_12, na.rm = T),
                            mean(cov50_int_13, na.rm = T),
                            mean(cov50_int_14, na.rm = T),
                            mean(cov50_int_15, na.rm = T),
                            mean(cov50_int_16, na.rm = T),
                            mean(cov50_int_17, na.rm = T),
                            mean(cov50_int_18, na.rm = T),
                            mean(cov50_int_19, na.rm = T),
                            mean(cov50_int_20, na.rm = T),
                            mean(cov50_int_21, na.rm = T),
                            mean(cov100_int_1, na.rm = T),
                            mean(cov100_int_2, na.rm = T),
                            mean(cov100_int_3, na.rm = T),
                            mean(cov100_int_4, na.rm = T),
                            mean(cov100_int_5, na.rm = T),
                            mean(cov100_int_6, na.rm = T),
                            mean(cov100_int_7, na.rm = T),
                            mean(cov100_int_8, na.rm = T),
                            mean(cov100_int_9, na.rm = T),
                            mean(cov100_int_10, na.rm = T),
                            mean(cov100_int_11, na.rm = T),
                            mean(cov100_int_12, na.rm = T),
                            mean(cov100_int_13, na.rm = T),
                            mean(cov100_int_14, na.rm = T),
                            mean(cov100_int_15, na.rm = T),
                            mean(cov100_int_16, na.rm = T),
                            mean(cov100_int_17, na.rm = T),
                            mean(cov100_int_18, na.rm = T),
                            mean(cov100_int_19, na.rm = T),
                            mean(cov100_int_20, na.rm = T),
                            mean(cov100_int_21, na.rm = T),
                            rmse20_int, rmse50_int, rmse100_int),
               "Metric" = rep(c("Bias", "Coverage", "RMSE"), each = 63),
               "Hline" = rep(c(0, 0.95, 0), each = 63),
               "Method" = factor(rep(c("CCA",
                                          rep(c("SI", "MI", "MMI", "B-MMI"),
                                       each = 5)), 9),
                                    levels = c("CCA", "SI", "MI", 
                                               "MMI", "B-MMI")),
               "Model" = rep(c("NA", rep(c("X+A+Y", "X+A*Y", "X*A+Y",
                                           "X*A+Y*A", "X*A*Y"),
                                           4)), 9),
               "Num_Clusters" = rep(rep(c(20, 50, 100), each = 21), 3))
  
  return(fig_dat_int)
  
}

figdat1 <- makeSimFig("out", 1000, 0.1, 0.1)
ggplot(data = figdat1[figdat1$Metric == "Bias", ],
       aes(x = Num_Clusters, y = Output, color = Model,
                           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Method, nrow = 1, scales = "fixed") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("Average bias") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines"))

ggplot(data = figdat1[figdat1$Metric == "Coverage", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Method, nrow = 1, scales = "fixed") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("Average coverage") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines"))

ggplot(data = figdat1[figdat1$Metric == "RMSE", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Method, nrow = 1, scales = "fixed") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("Average RMSE") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines"))



figdat2 <- makeSimFig("out2", 1000, 0.1, 0.1)
ggplot(data = figdat2[figdat2$Metric == "Bias", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Method, nrow = 1, scales = "fixed") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("Average bias") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines"))

ggplot(data = figdat2[figdat2$Metric == "Coverage", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Method, nrow = 1, scales = "fixed") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("Average coverage") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines"))

ggplot(data = figdat2[figdat2$Metric == "RMSE", ],
       aes(x = Num_Clusters, y = Output, color = Model,
           shape = Model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Method, nrow = 1, scales = "fixed") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  ylab("Average RMSE") +
  theme_light() +
  theme(panel.spacing = unit(0.75, "lines"))



# Old figure code
makeSimFig <- function(type, nsims, ICC_out, ICC_mod) {
  
  # 20 clusters
  num_clusters <- 20
  results <- array(NA, dim = c(nsims, 25))
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
  
  bias20 <- colMeans(results[, 1:9], na.rm = T) + 0.75
  
  cov20_1 <- 1*(results[, 1] - 1.96*results[, 10] < true_int &
                  results[, 1] + 1.96*results[, 10] > true_int)
  cov20_2 <- 1*(results[, 2] - 1.96*results[, 11] < true_int &
                  results[, 2] + 1.96*results[, 11] > true_int)
  cov20_3 <- 1*(results[, 3] - results[, 19]*results[, 12] < true_int &
                  results[, 3] + results[, 19]*results[, 12] > true_int)
  cov20_4 <- 1*(results[, 4] - results[, 20]*results[, 13] < true_int &
                  results[, 4] + results[, 20]*results[, 13] > true_int)
  cov20_5 <- 1*(results[, 5] - results[, 21]*results[, 14] < true_int &
                  results[, 5] + results[, 21]*results[, 14] > true_int)
  cov20_6 <- 1*(results[, 6] - results[, 22]*results[, 15] < true_int &
                  results[, 6] + results[, 22]*results[, 15] > true_int)
  cov20_7 <- 1*(results[, 7] - results[, 23]*results[, 16] < true_int &
                  results[, 7] + results[, 23]*results[, 16] > true_int)
  cov20_8 <- 1*(results[, 8] - results[, 24]*results[, 17] < true_int &
                  results[, 8] + results[, 24]*results[, 17] > true_int)
  cov20_9 <- 1*(results[, 9] - results[, 25]*results[, 18] < true_int &
                  results[, 9] + results[, 25]*results[, 18] > true_int)
  
  # Now with 50 clusters
  num_clusters <- 50
  results <- array(NA, dim = c(nsims, 25))
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
  
  bias50 <- colMeans(results[, 1:9], na.rm = T) + 0.75

  cov50_1 <- 1*(results[, 1] - 1.96*results[, 10] < true_int &
                  results[, 1] + 1.96*results[, 10] > true_int)
  cov50_2 <- 1*(results[, 2] - 1.96*results[, 11] < true_int &
                  results[, 2] + 1.96*results[, 11] > true_int)
  cov50_3 <- 1*(results[, 3] - results[, 19]*results[, 12] < true_int &
                  results[, 3] + results[, 19]*results[, 12] > true_int)
  cov50_4 <- 1*(results[, 4] - results[, 20]*results[, 13] < true_int &
                  results[, 4] + results[, 20]*results[, 13] > true_int)
  cov50_5 <- 1*(results[, 5] - results[, 21]*results[, 14] < true_int &
                  results[, 5] + results[, 21]*results[, 14] > true_int)
  cov50_6 <- 1*(results[, 6] - results[, 22]*results[, 15] < true_int &
                  results[, 6] + results[, 22]*results[, 15] > true_int)
  cov50_7 <- 1*(results[, 7] - results[, 23]*results[, 16] < true_int &
                  results[, 7] + results[, 23]*results[, 16] > true_int)
  cov50_8 <- 1*(results[, 8] - results[, 24]*results[, 17] < true_int &
                  results[, 8] + results[, 24]*results[, 17] > true_int)
  cov50_9 <- 1*(results[, 9] - results[, 25]*results[, 18] < true_int &
                  results[, 9] + results[, 25]*results[, 18] > true_int)

  # Now with 100 clusters
  num_clusters <- 100
  results <- array(NA, dim = c(nsims, 25))
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
  
  bias100 <- colMeans(results[, 1:9], na.rm = T) + 0.75

  cov100_1 <- 1*(results[, 1] - 1.96*results[, 10] < true_int &
                  results[, 1] + 1.96*results[, 10] > true_int)
  cov100_2 <- 1*(results[, 2] - 1.96*results[, 11] < true_int &
                  results[, 2] + 1.96*results[, 11] > true_int)
  cov100_3 <- 1*(results[, 3] - results[, 19]*results[, 12] < true_int &
                  results[, 3] + results[, 19]*results[, 12] > true_int)
  cov100_4 <- 1*(results[, 4] - results[, 20]*results[, 13] < true_int &
                  results[, 4] + results[, 20]*results[, 13] > true_int)
  cov100_5 <- 1*(results[, 5] - results[, 21]*results[, 14] < true_int &
                  results[, 5] + results[, 21]*results[, 14] > true_int)
  cov100_6 <- 1*(results[, 6] - results[, 22]*results[, 15] < true_int &
                  results[, 6] + results[, 22]*results[, 15] > true_int)
  cov100_7 <- 1*(results[, 7] - results[, 23]*results[, 16] < true_int &
                  results[, 7] + results[, 23]*results[, 16] > true_int)
  cov100_8 <- 1*(results[, 8] - results[, 24]*results[, 17] < true_int &
                  results[, 8] + results[, 24]*results[, 17] > true_int)
  cov100_9 <- 1*(results[, 9] - results[, 25]*results[, 18] < true_int &
                  results[, 9] + results[, 25]*results[, 18] > true_int)

  # Summary figure
  fig_dat <-
    data.frame("Output" = c(bias20[c(4:9)], bias50[c(4:9)],
                            bias100[c(4:9)],
                            mean(cov20_4, na.rm = T),
                            mean(cov20_5, na.rm = T),
                            mean(cov20_6, na.rm = T),
                            mean(cov20_7, na.rm = T),
                            mean(cov20_8, na.rm = T),
                            mean(cov20_9, na.rm = T),
                            mean(cov50_4, na.rm = T),
                            mean(cov50_5, na.rm = T),
                            mean(cov50_6, na.rm = T),
                            mean(cov50_7, na.rm = T),
                            mean(cov50_8, na.rm = T),
                            mean(cov50_9, na.rm = T),
                            mean(cov100_4, na.rm = T),
                            mean(cov100_5, na.rm = T),
                            mean(cov100_6, na.rm = T),
                            mean(cov100_7, na.rm = T),
                            mean(cov100_8, na.rm = T),
                            mean(cov100_9, na.rm = T)),
               "Metric" = rep(c("Bias", "Coverage"), each = 18),
               "Hline" = rep(c(0, 0.95), each = 18),
               "Method" = rep(c("X + Y + A", "X + A*Y", "X*A + Y",
                                "X*A + Y*A", "X*Y*A", "jomo"),
                              6),
               "Num_Clusters" = rep(rep(c(20, 50, 100), each = 6), 2))
  
  return(fig_dat)
  
}



figdat1 <- makeSimFig("out_XAM", 500, 0.1, 0.1)
ggplot(data = figdat1, aes(x = Num_Clusters, y = Output, color = Method,
                           shape = Method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Metric, nrow = 2, scales = "free") +
  geom_hline(aes(yintercept = Hline), color = "black") +
  #ggtitle("Modifier missing, ICC_out = ICC_mod = 0.1") +
  theme_light()
