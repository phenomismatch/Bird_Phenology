#!/bin/bash

#SBATCH --job-name=lc-br-Troglodytes_aedon
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=8G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_2019-01-16/lc-br-Troglodytes_aedon.out #STDOUT
#SBATCH -e /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_2019-01-16/lc-br-Troglodytes_aedon.err #STDERR

#echos name of node
echo `hostname`

#load singularity and run R script using singularity
module load singularity/3.0.2
singularity exec -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/8-logit-cubic-breeding/8-logit-cubic-breeding.R Troglodytes_aedon

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
