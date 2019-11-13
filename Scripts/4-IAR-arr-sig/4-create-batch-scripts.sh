#!/bin/bash

DATE="2019-11-13"

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  echo "#!/bin/bash

#SBATCH --job-name=sig-$temp
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 6 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=15G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/sig_output_$DATE/$temp-sig.out #STDOUT
#SBATCH -e /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/sig_output_$DATE/$temp-sig.err #STDERR

#echos name of node
echo \`hostname\`

#load R module and run script
#module load R/3.5.2
#Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/4-IAR-arr/4-IAR-arr.R $temp

#until singularity is sorted
module load singularity/3.0.2
singularity exec -B /labs/Tingley -B /UCHC /isg/shared/apps/R/3.5.2/R.sif Rscript /labs/Tingley/phenomismatch/Bird_Phenology/Scripts/4-IAR-arr-sig/4-IAR-arr-sig.R $temp

#displays amount of memory used
sstat --format=\"AveCPU,AvePages,AveRSS,MaxRSS,AveVMSize,MaxVMSize\" \$SLURM_JOBID.batch" > "species/4-$temp.sh"
done < ../../Data/test_species_list.txt
