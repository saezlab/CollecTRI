#!/bin/bash

#SBATCH --job-name=NTNU_knockTF
#SBATCH --output=NTNU_knockTF.out
#SBATCH --error=NTNU_knockTF.err
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=50G
#SBATCH -o NTNU_knockTF.out
#SBATCH -e NTNU_knockTF.err
#SBATCH -t 14:00:00

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleRabs
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/test_networks_knockTF_cluster_2sacct.R
