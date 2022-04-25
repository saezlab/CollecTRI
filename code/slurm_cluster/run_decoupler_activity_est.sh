#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --mem=15000
#SBATCH --job-name="decouplerAct"
#SBATCH --output=decouplerAct_%A_%a.out
#SBATCH --error=decouplerAct_%A_%a.err
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --signal=2
#SBATCH --array=1-14%7
#SBATCH --export=ALL
#SBATCH --mail-user=sophia.mueller-dott@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue

source ~/.bashrc
conda activate /net/data.isilon/ag-saez/bq_smueller/SOFTWARE/miniconda3/envs/decoupleR
Rscript /net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR/code/decoupler_activity_est.R $SLURM_ARRAY_TASK_ID
