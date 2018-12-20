#!/bin/bash

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
done < ../../../Data/eBird_species_list.txt
#done < ../../../Data/test_species_list.txt

DATE=`date +%Y-%m-%d`
cp 2-logit-cubic.R ../../../Data/Processed/halfmax_species/2-logit-cubic-$DATE.R
