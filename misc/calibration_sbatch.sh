#!/bin/bash
#SBATCH --time=00:25:00
#SBATCH --account=def-bolker
#SBATCH --job-name=ont_cal_NelderMead
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400M
module load r/4.1.0

loc=$(pwd)
ssh gra-login1 "cd $loc && git add --all && git commit -am 'calibration sync [skip ci]' && git pull"
#touch test_file.txt
Rscript ont_cal_NelderMead.R							
ssh gra-login1 "cd $loc && git pull && bash git_push.sh"
