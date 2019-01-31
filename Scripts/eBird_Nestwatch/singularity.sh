#!/bin/bash

grep libseccomp.so.2 * -lR > tfile

while read name
do
  temp="${name%.*}"
  cat ${temp}.out
done < tfile

rm tfile
