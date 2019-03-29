#!/bin/bash

#SBATCH --job-name=hm-arr-Setophaga_virens
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=25G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_species_2019-03-28/hm-arr-Setophaga_virens.out #STDOUT
#SBATCH -e /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_species_2019-03-28/hm-arr-Setophaga_virens.err #STDERR

#echos name of node
echo `hostname`

#load singularity and run R script using singularity
module load singularity/3.0.2
singularity exec -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/2-halfmax-arr/2-halfmax-arr.R Setophaga_virens

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
