#!/bin/bash
#BSUB -J R_job[1-500]
#BSUB -eo R-%I.e
module load R/4.0.2
Rscript crt-hte-missing-modifier.R