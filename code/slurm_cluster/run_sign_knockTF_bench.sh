#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=30:00:00
#SBATCH --mem=100000
#SBATCH --job-name="knockTF_test_sign"
#SBATCH --output=knockTF_test_sign.out
#SBATCH --mail-user=sophia.mueller-dott@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --cpus-per-task 64

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/test_sign_knockTF_cluster.R
