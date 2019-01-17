#!/bin/bash

#SBATCH --job-name=iar-Passerina_caerulea
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=8G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_2019-01-16/iar-Passerina_caerulea.out #STDOUT
#SBATCH -e /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_2019-01-16/iar-Passerina_caerulea.err #STDERR

#echos name of node
echo `hostname`

#load singularity module
#module load singularity/3.0.2
module load R/3.5.2

#run R script using singularity - species arg
#singularity exec /home/CAM/cyoungflesh/R.sif Rscript /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/4-IAR/4-IAR-model.R Passerina_caerulea
Rscript /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/4-IAR/4-IAR-model.R Passerina_caerulea

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
