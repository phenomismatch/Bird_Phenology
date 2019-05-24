#!/bin/bash

for a in *.txt
do
cat $a | head -n 12 | tail -n 2
done