#!/bin/bash

#argument is diretory of where to search for singularity error message
#script finds all STDERR that has the singularity error +
#print which node they occurred on +
#creates a text file with lines to resbumit jobs +
#removes all STDOUT SRDERR files associated with those jobs

grep libseccomp.so.2 $1 -lR > tfile

while read name
do
  temp="${name%.*}"
  head -n 1 ${temp}.out
  NAME2="${temp##*/}"
  NAME3="${NAME2##*-}"
  NAME4="${NAME3%.*}" 
  echo "sbatch 7-$NAME4.sh" >> miss_jobs.txt
  rm $temp*
done < tfile

num=$(wc -l tfile)
echo  "$num jobs have this issue"
rm tfile
