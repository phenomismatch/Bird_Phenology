#!/bin/bash

#SBATCH --job-name=iar-Petrochelidon_pyrrhonota
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=8G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_2019-01-16/iar-Petrochelidon_pyrrhonota.out #STDOUT
#SBATCH -e /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_2019-01-16/iar-Petrochelidon_pyrrhonota.err #STDERR

#echos name of node
echo `hostname`

#load R module (auto loads singularity and run script)
#module load R/3.5.2
#Rscript /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/4-IAR/4-IAR-model.R Petrochelidon_pyrrhonota

#load singularity and run R script using singularity
module load singularity/3.0.2
#singularity exec /home/CAM/cyoungflesh/R.sif Rscript /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/4-IAR/4-IAR-model.R Petrochelidon_pyrrhonota
singularity exec -B /UCHC /isg/shared/apps/R/3.5.2/R.sif /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/4-IAR/4-IAR-model.R Petrochelidon_pyrrhonota


#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
