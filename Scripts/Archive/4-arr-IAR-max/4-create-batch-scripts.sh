#!/bin/bash

DATE="2020-06-08"

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  echo "#!/bin/bash

#SBATCH --job-name=max-$DATE-$temp
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=15G #memory requested
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/arrival_IAR_max_$DATE/$temp-iar-max.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/arrival_IAR_max_$DATE/$temp-iar-max.err #STDERR

#echos name of node
echo \`hostname\`

module load gcc/6.4.0
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/4-arr-IAR-max/4-arr-IAR-max.R $temp 5000

#displays amount of memory used
sstat --format=\"AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize\" \$SLURM_JOBID.batch" > "species/4-$temp.sh"
done < ../../Data/IAR_species_list.txt
