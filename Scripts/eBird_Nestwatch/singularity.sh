#!/bin/bash

grep libseccomp.so.2 /home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_2019-01-30 -lR > tfile

while read name
do
  temp="${name%.*}"
  echo "$temp"
  head -n 1 ${temp}.out
done < tfile

num=$(wc -l tfile)
echo  "$num jobs have this issue"
rm tfile
