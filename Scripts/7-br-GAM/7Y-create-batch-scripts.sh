#!/bin/bash

DATE="2021-03-29"

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  echo "#!/bin/bash

#$ -cwd                         #directory (current working dir)
#$ -o br-gam-$temp-Y.joblog               #jobname
#$ -j y                         #combine STDOUT STDERR
#$ -pe shared 4                 #number of cores - for entire node: `-l exclusive; -pe node 1`
#$ -M cyoungl@mail              #mail address
#$ -m ea                        #email at end and abort times
#$ -l h_data=8G,h_rt=24:00:00  #resource request - run time in hours

#load modules
source /u/local/Modules/default/init/modules.sh
module load R/3.6.1

echo \`hostname\`

#export lib
export R_LIBS=/u/home/c/cyoungfl/R/x86_64-pc-linux-gnu-library/3.6

#run script - time call return info on memory usage
/usr/bin/time -apv Rscript /u/home/c/cyoungfl/Bird_Phenology/Scripts/7-br-GAM/7-br-GAM.R $temp Y

" > "species/a7-$temp-Y.sh"
done < ../../Data/arr_species_list.txt
