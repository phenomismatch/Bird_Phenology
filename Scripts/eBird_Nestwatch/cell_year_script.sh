#!/bin/bash

for a in *.out
do
    cat $a | tail -n 400 | grep species | tail -n 1
    JN=ls $a | head -c20
    qs | grep $JN
done
