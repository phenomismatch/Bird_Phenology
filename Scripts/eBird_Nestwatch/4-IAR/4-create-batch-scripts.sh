#!/bin/bash

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  DATE="2019-01-16"
  echo "#!/bin/bash

#SBATCH --job-name=iar-$temp
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=8G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_$DATE/iar-$temp.out #STDOUT
#SBATCH -e /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_$DATE/iar-$temp.err #STDERR

#echos name of node
echo \`hostname\`

#load R module (auto loads singularity and run script)
#module load R/3.5.2
#Rscript /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/4-IAR/4-IAR-model.R $temp

#load singularity and run R script using singularity
module load singularity/3.0.2
singularity exec /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Scripts/eBird_Nestwatch/4-IAR/4-IAR-model.R $temp

#displays amount of memory used
sstat --format=\"AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize\" \$SLURM_JOBID.batch" > "species/4-$temp.sh"
done < ../../../Data/IAR_species_list.txt
