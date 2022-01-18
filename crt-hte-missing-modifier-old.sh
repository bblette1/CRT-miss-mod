#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g

srun R CMD BATCH --vanilla "--args seed$SLURM_ARRAY_TASK_ID" crt-hte-missing-modifier.R crt-hte-missing-modifier-$SLURM_ARRAY_TASK_ID.Rout