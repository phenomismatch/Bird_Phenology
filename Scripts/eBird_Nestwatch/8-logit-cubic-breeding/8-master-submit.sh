#!/bin/bash

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/8-$temp.sh
done < ../../../Data/IAR_species_list.txt
#done < ../../../Data/test_species_list.txt

DATE="2019-01-16"
cp 8-logit-cubic-breeding.R ../../../Data/Processed/halfmax_breeding_$DATE/8-logit-cubic-breeding-$DATE.R
