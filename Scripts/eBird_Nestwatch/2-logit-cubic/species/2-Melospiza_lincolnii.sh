#!/bin/bash

#SBATCH --job-name=lc-arr-Melospiza_lincolnii
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=8G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_species_2019-01-16/lc-arr-Melospiza_lincolnii.out #STDOUT
#SBATCH -e /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_species_2019-01-16/lc-arr-Melospiza_lincolnii.err #STDERR

#echos name of node
echo `hostname`

#load R module (auto loads singularity and run script)
#module load R/3.5.2
#Rscript /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/2-logit-cubic/2-logit-cubic.R Melospiza_lincolnii

#load singularity and run R script using singularity
module load singularity/3.0.2
singularity exec /home/CAM/cyoungflesh/R.sif Rscript /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/2-logit-cubic/2-logit-cubic.R Melospiza_lincolnii

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
