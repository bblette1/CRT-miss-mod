# CRT-miss-mod

This repository contains files used for the simulations and applications of "Assessing treatment effect heterogeneity in the presence of missing effect modifier data in cluster-randomized trials" (submitted manuscript).

### Replication of simulation study

The file *crt-hte-missing-modifier.R* can be used to replicate Scenario 1 of the simulation study and the file *crt-hte-missing-modifier.R* can be used to replicate Scenario 2 of the simulation study. Run the simulator function in each file for 2000 replications for each of `num_clusters %in% c(20, 50, 100)` with the input seed structure described at the bottom of the files. To make the main figures of the paper, see *make-sim-fig.R*.

### Replication of application

To replicate the data analysis presented in the paper, you must download publically available data from the Work, Family, and Health Study [https://www.icpsr.umich.edu/web/DSDR/studies/36158/versions/V2](https://www.icpsr.umich.edu/web/DSDR/studies/36158/versions/V2). Then use the code file *data-analysis.R*.
