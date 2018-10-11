#!/bin/bash

#SBATCH --job-name=logit-cubic
#SBATCH -N 1 #number of tasks 
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --qos=general #queue (same as partition)
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mem=50G #memory requested
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o %j.out #STDOUT
#SBATCH -e %j.err #STDERR

echo `hostname` #echos name of node

module load R/3.5.1

export R_LIBS=/home/CAM/cyoungflesh/R_libs

Rscript 2-logit-cubic.R