#!/bin/bash
#SBATCH -J matlab
#SBATCH -o matlab.out
#SBATCH -N 2
#SBATCH -c 8
#SBATCH --mem=64G
module load Matlab/R2020a
matlab -nodisplay -nodesktop -r "run /hpc/group/katzlab/POI/Analysis_v1/complete_analysis.m"