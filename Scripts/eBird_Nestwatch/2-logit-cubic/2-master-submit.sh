#!/bin/bash

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
#done < ../../../Data/eBird_species_list.txt
done < ../../../Data/test_species_list.txt
