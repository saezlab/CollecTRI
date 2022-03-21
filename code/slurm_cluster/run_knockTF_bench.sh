#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=25:00:00
#SBATCH --mem=15000
#SBATCH --job-name="knockTF_test_networks"
#SBATCH --output=knockTF_test_networks.out
#SBATCH --mail-user=sophia.mueller-dott@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --cpus-per-task 64

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/test_networks_knockTF_cluster.R
