#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=25:00:00
#SBATCH --mem=15000
#SBATCH --job-name="dorothea_test_networks"
#SBATCH --output=dorothea_test_networks.out
#SBATCH --mail-user=sophia.mueller-dott@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --cpus-per-task 64

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/test_networks_dorothea_cluster.R
