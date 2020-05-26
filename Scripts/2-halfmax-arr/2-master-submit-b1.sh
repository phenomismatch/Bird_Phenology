#!/bin/bash

DATE="2020-05-07"

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
done < ../../Data/eBird_species_list_b1.txt

cp 2-halfmax-arr.R ../../Data/Processed/halfmax_arrival_$DATE/2-halfmax-arr-$DATE.R
