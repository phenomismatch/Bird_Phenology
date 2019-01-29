#!/bin/bash

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../../Data/IAR_species_list.txt
#done < ../../../Data/test_species_list.txt

DATE="2019-01-16"
cp 4-IAR-model.R ../../../Data/Processed/IAR_output_$DATE/4-IAR-model-$DATE.R
