#!/bin/bash

DATE="2020-02-26"

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
done < ../../Data/eBird_species_list_b2.txt
