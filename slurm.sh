#!/bin/bash
#SBATCH --job-name=kura64
#SBATCH --workdir=/master/home/qiqi/git/Diamond1D
#SBATCH --output=kura64.out
#SBATCH --error=kura64.err
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
 
make run
