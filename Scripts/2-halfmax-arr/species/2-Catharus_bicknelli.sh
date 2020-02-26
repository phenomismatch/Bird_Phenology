#!/bin/bash

#SBATCH --job-name=hm-arr-Catharus_bicknelli
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=25G #memory requested
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_arrival_2020-02-26/hm-arr-Catharus_bicknelli.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_arrival_2020-02-26/hm-arr-Catharus_bicknelli.err #STDERR

#echos name of node
echo `hostname`

#load R module and run script
#module load R/3.5.2
#Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/2-halfmax-arr/2-halfmax-arr.R Catharus_bicknelli

#until singularity is sorted
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/2-halfmax-arr/2-halfmax-arr.R Catharus_bicknelli

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
