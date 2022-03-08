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
