#!/bin/bash

#SBATCH --job-name=knockTF
#SBATCH -t 20:00:00
#SBATCH -n 4
#SBATCH -p multi
#SBATCH -N 4
#SBATCH --output knockTF.out
#SBATCH --error knockTF.err

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/test_networks_knockTF_cluster.R
