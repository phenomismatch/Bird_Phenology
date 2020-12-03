#!/bin/bash

#SBATCH --job-name=br-gam-Quiscalus_quiscula
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=25G #memory requested
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/breeding_GAM_2020-12-03/br-gam-Quiscalus_quiscula.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/breeding_GAM_2020-12-03/br-gam-Quiscalus_quiscula.err #STDERR

#echos name of node
echo `hostname`

module load gcc/6.4.0
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/8-br-GAM/8-br-GAM.R Quiscalus_quiscula

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
