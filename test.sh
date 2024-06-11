#!/bin/bash
#SBATCH -J doe
#SBATCH -o doe.out
#SBATCH -N 2
#SBATCH -c 8
#SBATCH --mem=64G
module load Matlab/R2022b
matlab -nodisplay -nodesktop -r "run /hpc/group/katzlab/POI/Analysis_v1/complete_analysis.m"