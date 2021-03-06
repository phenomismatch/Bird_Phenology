#!/bin/bash

#$ -cwd                         #directory (current working dir)
#$ -o br-gam-Tyrannus_tyrannus-E.joblog               #jobname
#$ -j y                         #combine STDOUT STDERR
#$ -pe shared 4                 #number of cores - for entire node: 
#$ -M cyoungl@mail              #mail address
#$ -m ea                        #email at end and abort times
#$ -l h_data=8G,h_rt=24:00:00  #resource request - run time in hours

#load modules
source /u/local/Modules/default/init/modules.sh
module load R/3.6.1

echo `hostname`

#export lib
export R_LIBS=/u/home/c/cyoungfl/R/x86_64-pc-linux-gnu-library/3.6

#run script - time call return info on memory usage
/usr/bin/time -apv Rscript /u/home/c/cyoungfl/Bird_Phenology/Scripts/7-br-GAM/7-br-GAM.R Tyrannus_tyrannus E


