#!/bin/bash
#SBATCH --time=00:01:00
#SBATCH --account=def-bolker
#SBATCH --job-name=git_test
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=128M


ssh gra-login1 < bash git_push.sh $(pwd)