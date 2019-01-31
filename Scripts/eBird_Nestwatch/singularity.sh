#!/bin/bash

grep libseccomp.so.2 $1 -lR > tfile

while read name
do
  temp="${name%.*}"
  head -n 1 ${temp}.out
done < tfile

num=$(wc -l tfile)
echo  "$num jobs have this issue"
rm tfile
