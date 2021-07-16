#!/bin/bash
#SBATCH --time=00:01:00
#SBATCH --account=def-bolker
#SBATCH --job-name=ont_cal_NelderMead
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400M
module load r/4.1.0
#00:02:00
# Rscript ont_cal_NelderMead.R
loc=$(pwd)
ssh gra-login1 < git_push.sh $loc
