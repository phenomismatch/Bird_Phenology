#!/bin/bash

#SBATCH --job-name=arr-GAM-Hirundo_rustica
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=40G #memory requested
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/arrival_GAM_2020-07-10/arr-GAM-Hirundo_rustica.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/arrival_GAM_2020-07-10/arr-GAM-Hirundo_rustica.err #STDERR

#echos name of node
echo `hostname`

module load gcc/6.4.0
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/2-arr-GAM/2-arr-GAM.R Hirundo_rustica 2002 2017

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
