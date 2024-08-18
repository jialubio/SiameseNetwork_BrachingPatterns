#!/bin/bash
#SBATCH --output=matlab.out
#SBATCH --array=1-100
#SBATCH -p scavenger	
#SBATCH --mem=2G
 /opt/apps/rhel7/matlabR2017b/bin/matlab -nodesktop -nodisplay -singleCompThread -r "rank=$SLURM_ARRAY_TASK_ID;OptimalPattern_homogeneous;quit"
