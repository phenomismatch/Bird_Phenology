#!/bin/bash

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
#done < ../../../Data/eBird_species_list.txt
done < ../../../Data/test_species_list.txt

DATE="2019-01-16"
cp 2-logit-cubic.R ../../../Data/Processed/halfmax_species_$DATE/2-logit-cubic-$DATE.R
