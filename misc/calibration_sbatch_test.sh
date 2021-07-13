#!/bin/bash
#SBATCH --time=00:00:15
#SBATCH --account=def-bolker
#SBATCH --job-name=git_test
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64M

touch calibration_sbatch_test.txt 
loc=$(pwd)
ssh gra-login2 "cd $loc; echo $(pwd); bash git_push.sh"