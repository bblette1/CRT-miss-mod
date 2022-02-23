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
  se20 <- colMeans(results[, 10:18], na.rm = T)
  
  cov20_1 <- 1*(results[, 1] - 1.96*results[, 10] < -0.75 &
                  results[, 1] + 1.96*results[, 10] > -0.75)
  cov20_2 <- 1*(results[, 2] - 1.96*results[, 11] < -0.75 &
                  results[, 2] + 1.96*results[, 11] > -0.75)
  cov20_3 <- 1*(results[, 3] - results[, 19]*results[, 12] < -0.75 &
                  results[, 3] + results[, 19]*results[, 12] > -0.75)
  cov20_4 <- 1*(results[, 4] - results[, 20]*results[, 13] < -0.75 &
                  results[, 4] + results[, 20]*results[, 13] > -0.75)
  cov20_5 <- 1*(results[, 5] - results[, 21]*results[, 14] < -0.75 &
                  results[, 5] + results[, 21]*results[, 14] > -0.75)
  cov20_6 <- 1*(results[, 6] - results[, 22]*results[, 15] < -0.75 &
                  results[, 6] + results[, 22]*results[, 15] > -0.75)
  cov20_7 <- 1*(results[, 7] - results[, 23]*results[, 16] < -0.75 &
                  results[, 7] + results[, 23]*results[, 16] > -0.75)
  cov20_8 <- 1*(results[, 8] - results[, 24]*results[, 17] < -0.75 &
                  results[, 8] + results[, 24]*results[, 17] > -0.75)
  cov20_9 <- 1*(results[, 9] - results[, 25]*results[, 18] < -0.75 &
                  results[, 9] + results[, 25]*results[, 18] > -0.75)
  
  ser20_1 <- sd(results[, 1], na.rm = T) / mean(results[, 10], na.rm = T)
  ser20_2 <- sd(results[, 2], na.rm = T) / mean(results[, 11], na.rm = T)
  ser20_3 <- sd(results[, 3], na.rm = T) / mean(results[, 12], na.rm = T)
  ser20_4 <- sd(results[, 4], na.rm = T) / mean(results[, 13], na.rm = T)
  ser20_5 <- sd(results[, 5], na.rm = T) / mean(results[, 14], na.rm = T)
  ser20_6 <- sd(results[, 6], na.rm = T) / mean(results[, 15], na.rm = T)
  ser20_7 <- sd(results[, 7], na.rm = T) / mean(results[, 16], na.rm = T)
  ser20_8 <- sd(results[, 8], na.rm = T) / mean(results[, 17], na.rm = T)
  ser20_9 <- sd(results[, 9], na.rm = T) / mean(results[, 18], na.rm = T)
  
  
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
  se50 <- colMeans(results[, 10:18], na.rm = T)
  
  cov50_1 <- 1*(results[, 1] - 1.96*results[, 10] < -0.75 &
                  results[, 1] + 1.96*results[, 10] > -0.75)
  cov50_2 <- 1*(results[, 2] - 1.96*results[, 11] < -0.75 &
                  results[, 2] + 1.96*results[, 11] > -0.75)
  cov50_3 <- 1*(results[, 3] - results[, 19]*results[, 12] < -0.75 &
                  results[, 3] + results[, 19]*results[, 12] > -0.75)
  cov50_4 <- 1*(results[, 4] - results[, 20]*results[, 13] < -0.75 &
                  results[, 4] + results[, 20]*results[, 13] > -0.75)
  cov50_5 <- 1*(results[, 5] - results[, 21]*results[, 14] < -0.75 &
                  results[, 5] + results[, 21]*results[, 14] > -0.75)
  cov50_6 <- 1*(results[, 6] - results[, 22]*results[, 15] < -0.75 &
                  results[, 6] + results[, 22]*results[, 15] > -0.75)
  cov50_7 <- 1*(results[, 7] - results[, 23]*results[, 16] < -0.75 &
                  results[, 7] + results[, 23]*results[, 16] > -0.75)
  cov50_8 <- 1*(results[, 8] - results[, 24]*results[, 17] < -0.75 &
                  results[, 8] + results[, 24]*results[, 17] > -0.75)
  cov50_9 <- 1*(results[, 9] - results[, 25]*results[, 18] < -0.75 &
                  results[, 9] + results[, 25]*results[, 18] > -0.75)
  
  ser50_1 <- sd(results[, 1], na.rm = T) / mean(results[, 10], na.rm = T)
  ser50_2 <- sd(results[, 2], na.rm = T) / mean(results[, 11], na.rm = T)
  ser50_3 <- sd(results[, 3], na.rm = T) / mean(results[, 12], na.rm = T)
  ser50_4 <- sd(results[, 4], na.rm = T) / mean(results[, 13], na.rm = T)
  ser50_5 <- sd(results[, 5], na.rm = T) / mean(results[, 14], na.rm = T)
  ser50_6 <- sd(results[, 6], na.rm = T) / mean(results[, 15], na.rm = T)
  ser50_7 <- sd(results[, 7], na.rm = T) / mean(results[, 16], na.rm = T)
  ser50_8 <- sd(results[, 8], na.rm = T) / mean(results[, 17], na.rm = T)
  ser50_9 <- sd(results[, 9], na.rm = T) / mean(results[, 18], na.rm = T)
  
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
  se100 <- colMeans(results[, 10:18], na.rm = T)
  
  cov100_1 <- 1*(results[, 1] - 1.96*results[, 10] < -0.75 &
                  results[, 1] + 1.96*results[, 10] > -0.75)
  cov100_2 <- 1*(results[, 2] - 1.96*results[, 11] < -0.75 &
                  results[, 2] + 1.96*results[, 11] > -0.75)
  cov100_3 <- 1*(results[, 3] - results[, 19]*results[, 12] < -0.75 &
                  results[, 3] + results[, 19]*results[, 12] > -0.75)
  cov100_4 <- 1*(results[, 4] - results[, 20]*results[, 13] < -0.75 &
                  results[, 4] + results[, 20]*results[, 13] > -0.75)
  cov100_5 <- 1*(results[, 5] - results[, 21]*results[, 14] < -0.75 &
                  results[, 5] + results[, 21]*results[, 14] > -0.75)
  cov100_6 <- 1*(results[, 6] - results[, 22]*results[, 15] < -0.75 &
                  results[, 6] + results[, 22]*results[, 15] > -0.75)
  cov100_7 <- 1*(results[, 7] - results[, 23]*results[, 16] < -0.75 &
                  results[, 7] + results[, 23]*results[, 16] > -0.75)
  cov100_8 <- 1*(results[, 8] - results[, 24]*results[, 17] < -0.75 &
                  results[, 8] + results[, 24]*results[, 17] > -0.75)
  cov100_9 <- 1*(results[, 9] - results[, 25]*results[, 18] < -0.75 &
                  results[, 9] + results[, 25]*results[, 18] > -0.75)
  
  ser100_1<- sd(results[, 1], na.rm = T) / mean(results[, 10], na.rm = T)
  ser100_2<- sd(results[, 2], na.rm = T) / mean(results[, 11], na.rm = T)
  ser100_3<- sd(results[, 3], na.rm = T) / mean(results[, 12], na.rm = T)
  ser100_4<- sd(results[, 4], na.rm = T) / mean(results[, 13], na.rm = T)
  ser100_5<- sd(results[, 5], na.rm = T) / mean(results[, 14], na.rm = T)
  ser100_6<- sd(results[, 6], na.rm = T) / mean(results[, 15], na.rm = T)
  ser100_7<- sd(results[, 7], na.rm = T) / mean(results[, 16], na.rm = T)
  ser100_8<- sd(results[, 8], na.rm = T) / mean(results[, 17], na.rm = T)
  ser100_9<- sd(results[, 9], na.rm = T) / mean(results[, 18], na.rm = T)
  
  # Summary figure
  fig_dat <-
    data.frame("Output" = c(bias20[c(4:9)], bias50[c(4:9)],
                            bias100[c(4:9)], se20[c(4:9)],
                            se50[c(4:9)], se100[c(4:9)],
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
                            mean(cov100_9, na.rm = T),
                            ser20_4, ser20_5, ser20_6, ser20_7, ser20_8,
                            ser20_9, ser50_4, ser50_5, ser50_6, ser50_7,
                            ser50_8, ser50_9, ser100_4, ser100_5,
                            ser100_6, ser100_7, ser100_8, ser100_9),
               "Metric" = rep(c(" Bias", "ASE", "Coverage",
                                "ESE / ASE"), each = 18),
               "Hline" = rep(c(0, 0, 0.95, 1), each = 18),
               "Method" = rep(c("X + Y + A", "X + A*Y", "X*A + Y",
                                "X*A + Y*A", "X*Y*A", "BART"), 12),
               "Num_Clusters" = rep(rep(c(20, 50, 100), each = 6), 4))
  
  return(fig_dat)
  
}