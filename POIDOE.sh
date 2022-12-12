#!/bin/bash
#SBATCH -J doe
#SBATCH -o doe.out
#SBATCH --partition=katzlab
#SBATCH -n 92
#92

module load Matlab/R2022b

matlab -nodisplay -nodesktop -r "run /hpc/group/katzlab/POI/Analysis_v1/complete_analysis.m"