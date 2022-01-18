#!/bin/bash
#BSUB -J R_job[1-10]
#BSUB -oo R-%I.o
#BSUB -eo R-%I.e
module load R/3.5.0
Rscript crt-hte-missing-modifier.R