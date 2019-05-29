#!/bin/bash

for a in *.txt
do
cat $a | head -n 2
cat $a | head -n 12 | tail -n 2
cat $a | head -n 6 | tail -n 1
done