#!/bin/bash
#SBATCH --time=00:25:00
#SBATCH --account=def-bolker
#SBATCH --job-name=ont_cal_NelderMead
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400M
module load r/4.1.0
git pull
Rscript ont_cal_NelderMead.R

git config --global user.name 'calibration-bot'
git config --global user.email 'calibration-bot@users.noreply.github.com'
git add --all
git commit -am "Weekly regression testing calibration [skip ci]"
git push
