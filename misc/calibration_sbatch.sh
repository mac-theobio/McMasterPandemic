#!/bin/bash
#SBATCH --time=00:03:00
#SBATCH --account=def-bolker
#SBATCH --job-name=ont_cal_NelderMead
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400M
module load r/4.1.0
#00:02:00
# Rscript ont_cal_NelderMead.R
loc=$(pwd)
ssh gra-login1 cd $loc && git add --all && git commit -am "test"  && git pull
ssh gra-login2 cd $loc && bash git_push.sh 
