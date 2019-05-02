#!/bin/bash

for a in *.o*
do
cat $a | tail -n 400 | grep species | tail -n 1
JN=`ls $a | head -c15 | tail -c+3`
qstat -u cyoungflesh | grep $JN | tail -c8
done
