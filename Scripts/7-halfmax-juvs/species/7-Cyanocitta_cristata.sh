#!/bin/bash

#SBATCH --job-name=hm-juvs-Cyanocitta_cristata
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=25G #memory requested
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_juvs_2020-03-25/hm-juvs-Cyanocitta_cristata.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_juvs_2020-03-25/hm-juvs-Cyanocitta_cristata.err #STDERR

#echos name of node
echo `hostname`

#load R module and run script
# module load R/3.5.2
# Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/7-halfmax-juvs/7-halfmax-juvs.R Cyanocitta_cristata

#until singularity is sorted
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/7-halfmax-juvs/7-halfmax-juvs.R Cyanocitta_cristata

#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
