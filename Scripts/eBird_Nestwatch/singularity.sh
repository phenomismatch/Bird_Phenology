#!/bin/bash

grep libseccomp.so.2 /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_2019-01-30 -lR > tfile

while read name
do
  temp="${name%.*}"
  cat ${temp}.out
done < tfile

rm tfile
