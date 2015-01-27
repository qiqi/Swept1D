#!/bin/bash
#SBATCH --job-name=kuramotoDiamondTest
#SBATCH --workdir=/master/home/qiqi/git/cloud_pde/1d
#SBATCH --output=kuramotoDiamondTest4.out
#SBATCH --error=kuramotoDiamondTest4.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
 
make run
