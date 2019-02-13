#!/bin/bash

DATE="2019-02-14"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../../Data/IAR_species_list-2019-02-02.txt
#done < ../../../Data/test_species_list.txt

cp 4-IAR-model.R ../../../Data/Processed/IAR_output_$DATE/4-IAR-model-$DATE.R
