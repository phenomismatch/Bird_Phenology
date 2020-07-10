#!/bin/bash

for a in *.txt
do
cat $a | head -n 1
cat $a | head -n 16 | tail -n 3
printf '\n'
done