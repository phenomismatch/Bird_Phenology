#!/bin/bash

DATE="2021-03-29"
mkdir /u/home/c/cyoungfl/Bird_Phenology/Data/Processed/breeding_GAM_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  qsub species/a7-$temp-F.sh
done < ../../Data/arr_species_list.txt

cp 7-br-GAM.R ../../Data/Processed/breeding_GAM_$DATE/7-br-GAM-F-$DATE.R
