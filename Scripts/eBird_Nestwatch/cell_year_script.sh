#!/bin/bash

for a in *.out
do
    cat a | tail -n 400 | grep species | tail -n 1
done