#!/bin/bash

#SBATCH --job-name=lc-Icterus_galbula
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=10G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o lc-Icterus_galbula.out #STDOUT
#SBATCH -e lc-Icterus_galbula.err #STDERR

#echos name of node
echo `hostname`

#load R module
module load R/3.5.1

#R library
export R_LIBS=/home/CAM/cyoungflesh/R_libs

#run R script - species arg
Rscript 2-logit-cubic.R Icterus_galbula

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
