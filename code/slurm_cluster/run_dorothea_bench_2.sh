#!/bin/bash

#SBATCH --job-name=NTNU_dorotheabench
#SBATCH --output=NTNU_dorotheabench.out
#SBATCH --error=NTNU_dorotheabench.err
#SBATCH --export=ALL
#SBATCH --nodes=5
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=50G
#SBATCH -o NTNU_dorotheabench.out
#SBATCH -e NTNU_dorotheabench.err
#SBATCH -t 14:00:00

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/test_networks_cluster_2.R
