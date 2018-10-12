#!/bin/bash

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch $temp.sh
done < ~/Google_Drive/R/Bird_Phenology/Data/eBird_species_list.txt