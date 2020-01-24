#!/bin/bash

DATE="2019-10-31"

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  echo "#!/bin/bash

#SBATCH --job-name=hm-br-$temp
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=25G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_$DATE/hm-br-$temp.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_$DATE/hm-br-$temp.err #STDERR

#echos name of node
echo \`hostname\`

#load singularity and run R script using singularity
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/8-halfmax-br/8-halfmax-br.R $temp

#displays amount of memory used
sstat --format=\"AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize\" \$SLURM_JOBID.batch" > "species/8-$temp.sh"
done < ../../Data/IAR_species_list.txt