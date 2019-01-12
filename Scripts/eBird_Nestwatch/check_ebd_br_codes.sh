#!/bin/bash

#SBATCH --job-name=ebd
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --qos=general #queue (same as partition)
#SBATCH --partition=general #partition - can also specify 'himem'
#SBATCH --mem=10G #memory requested
#SBATCH --mail-type=END #when to send email (on job completion)
#SBATCH --mail-user=casey.youngflesh@uconn.edu #email address for notification
#SBATCH -o /home/CAM/cyoungflesh/ebd.out #STDOUT
#SBATCH -e /home/CAM/cyoungflesh/ebd.err #STDERR

#echos name of node
echo `hostname`

cd /home/CAM/cyoungflesh/phenomismatch/useful_datasets/eBird/ebd_relFeb-2018


awk -F'\t' '$6=="Contopus virens" {print $0}' ebd_relFeb-2018.txt > Contopus_virens.txt

awk -F'\t' '$11=="C3" || $11=="C4" {print $0}' Contopus_virens.txt > C34_Contopus_virens.txt
