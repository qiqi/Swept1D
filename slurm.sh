#!/bin/bash
#SBATCH --job-name=swept32
#SBATCH --workdir=/master/home/qiqi/git/Diamond1D
#SBATCH --output=swept32.out
#SBATCH --error=swept32.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
 
make run_heat
make run_kuramoto
