#!/bin/bash

#SBATCH --job-name=decouplerAct
#SBATCH -t 20:00:00
#SBATCH -n 4
#SBATCH -p multi
#SBATCH -N 4
#SBATCH --output decouplerAct.out
#SBATCH --error decouplerAct.err

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/decoupler_activity_est.R
