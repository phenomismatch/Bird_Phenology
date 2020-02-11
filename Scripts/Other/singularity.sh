#!/bin/bash

#first argument is search term
#second argument is diretory of where to search for singularity error message
#script finds all STDERR that has the singularity error +
#print which node they occurred on +
#creates a text file with lines to resbumit jobs +
#removes all STDOUT SRDERR files associated with those jobs

#SUM of two lines below should == total number of jobs submitted
#grep 'completed' *.out -lR | wc -l #jobs completed
#squeue -u cyoungflesh | grep lc-br | wc -l #lc-br jobs running

grep $1 $2 -lR > tfile

while read name
do
  temp="${name%.*}"
  head -n 1 ${temp}.out
  NAME2="${temp##*/}"
  NAME3="${NAME2##*-}"
  NAME4="${NAME3%.*}" 
  echo "sbatch 2-$NAME4.sh" >> miss_jobs.txt
  rm $temp*
done < tfile

num=$(wc -l tfile)
echo  "$num jobs have this issue"
rm tfile
