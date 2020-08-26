#!/bin/bash

DATE="2020-07-21"

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../Data/IAR_species_list_p2.txt
