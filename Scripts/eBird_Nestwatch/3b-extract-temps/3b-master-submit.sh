#!/bin/bash

while read year
do
  sbatch years/3b-$year.sh
done < daymet_years.txt
