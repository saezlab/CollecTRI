#!/bin/bash

#SBATCH --job-name=knockTFsign
#SBATCH -t 20:00:00
#SBATCH -n 4
#SBATCH -p multi
#SBATCH -N 4
#SBATCH --output knockTF_sign.out
#SBATCH --error knockTF_sign.err

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/test_sign_knockTF_cluster.R
