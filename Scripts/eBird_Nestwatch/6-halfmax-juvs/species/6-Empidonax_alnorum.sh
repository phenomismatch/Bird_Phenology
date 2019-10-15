#!/bin/bash

#SBATCH --job-name=hm-juvs-Empidonax_alnorum
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=25G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_juvs_2019-10-15/hm-juvs-Empidonax_alnorum.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_juvs_2019-10-15/hm-juvs-Empidonax_alnorum.err #STDERR

#echos name of node
echo `hostname`

#load R module and run script
# module load R/3.5.2
# Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/6-halfmax-juvs/6-halfmax-juvs.R Empidonax_alnorum

#until singularity is sorted
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/6-halfmax-juvs/6-halfmax-juvs.R Empidonax_alnorum


#displays amount of memory used
sstat --format="AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize" $SLURM_JOBID.batch
