#!/bin/bash

grep libseccomp.so.2 $1 -lR > tfile

while read name
do
  temp="${name%.*}"
  head -n 1 ${temp}.out
  NAME2="${temp##*/}"
  NAME3="${NAME2##*-}"
  NAME4="${NAME3%.*}" 
  echo "sbatch 7-$NAME4.sh" >> miss_jobs.txt
done < tfile

num=$(wc -l tfile)
echo  "$num jobs have this issue"
rm tfile
