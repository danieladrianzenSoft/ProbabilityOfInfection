#!/bin/bash
#SBATCH -J doe
#SBATCH -o doe.out
#SBATCH --partition=scavenger
#SBATCH --mem=64G
module load Matlab/R2022b
matlab -nodisplay -nodesktop -r "run /hpc/group/katzlab/POI/Cluster/complete_analysis.m"