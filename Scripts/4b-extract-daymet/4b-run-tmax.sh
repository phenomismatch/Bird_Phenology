#!/bin/bash

#SBATCH --job-name=extract-tmax
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 2 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=250G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/daymet/extract-tmax.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/daymet/extract-tmax.err #STDERR

#echos name of node
echo `hostname`

#load R module and run script
module load R/3.5.2
Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/4b-extract-daymet/4b-extract-tmax.R

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
